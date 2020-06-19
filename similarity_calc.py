import math
import numpy
import sys

class Peak:
    def __init__(self,mass,intensity):
        self.mass = mass
        self.intensity = intensity
        self.sqrt_intensity = math.sqrt(self.intensity)
    
    def __repr__(self):
        return f"{self.__class__.__name__}({self.mass},{self.intensity})"

    def compare_peaks(self, peak_2, tolerance):
        """Compares to another peak and returns true if masses are within given tolerance"""
        if (self.mass >= peak_2.mass-tolerance and self.mass <= peak_2.mass+tolerance):
            return True
        else:
            return False

class Spectrum:
    def __init__(self):
        self.peaks = [] #list of peaks in spectrum
        self.feature_id=None
        
    def __repr__(self):
        return f"{self.__class__.__name__}({self.peaks})"

    def add_peak(self,mass, intensity):
        self.peaks.append(Peak(mass,intensity))
    
    def set_peak_id(self,feature_id):
        self.feature_id = feature_id

    def compare_precursor_mass(self,spectrum_2,threshold):
        if ((self.pep_mass-spectrum_2.pep_mass)>(-threshold) and (self.pep_mass-spectrum_2.pep_mass)<threshold):
            return True
        else:
            return False
    
    def euclidean_scale(self):
        """Scales the intensity of each Peak object using Euclidean norm"""
        
        #get list of intensities from spectrum
        intensities=[]
        for peak in self.peaks:
            intensities.append(peak.sqrt_intensity)
        
        #calculate norm
        norm=numpy.linalg.norm(intensities)

        #set each peak's scaled intensity
        for peak in self.peaks:
            peak.norm_scaled=(peak.sqrt_intensity)/norm



def mgf_reader(file_path):
    """Takes the path to a .mgf file. Returns a list of Spectrum objects, each containing the attributes contained in the
    file and a list of Peak objects"""
    file = open(file_path)
    #split file into list of lines
    lines = file.read().split('\n')

    #initiate list to hold all the Spectrum objects
    spectra_list=[]

    for line in lines:
        if 'BEGIN IONS' in line:
            spectrum_lines = [] #initiate new spectra list to contain lines of the next spectra data
            spectrum = Spectrum() #initiate new spectrum object
            continue
   
        if 'END IONS' in line:
            #at end of current spectra, parse the data
            
            
            for i in spectrum_lines:

                #extract parameter values and assign to attributes
                if "=" in i:
                    parameter,value = i.split("=")
                    if "feature_id".casefold() in parameter.casefold():
                        spectrum.set_peak_id(value)
                    elif "pepmass".casefold() in parameter.casefold():
                        spectrum.pep_mass = float(value)
                    elif "scans".casefold() in parameter.casefold():
                        spectrum.scans = float(value)
                    elif "RTinseconds".casefold() in parameter.casefold():
                        spectrum.retention_time = float(value)
                    elif "charge".casefold() in parameter.casefold():
                        spectrum.charge = value
                    elif "mslevel".casefold() in parameter.casefold():
                        spectrum.ms_level = int(value)

                #when it reaches a line with no '=' then it is one of the peaks, which is added to the spectrum
                else:
                    mass, intensity = i.split()
                    mass = float(mass)
                    intensity = float(intensity)
                    spectrum.add_peak(mass, intensity)

            #intensities scaled by Euclidean norm
            spectrum.euclidean_scale()

            spectra_list.append(spectrum)  
         
        #add line to existing spectra info
        spectrum_lines.append(line)

    file.close()
    return spectra_list


def cosine_score(spectrum_one, spectrum_two, tolerance):
    """takes two MS2 Spectrum objects and returns the cosine similarity score"""  
    peak_pairs=[]
    #make list of pairs of matched peaks within given tolerance
    for peak1 in spectrum_one.peaks:
        for peak2 in spectrum_two.peaks:
            if peak1.compare_peaks(peak2,tolerance):
                product=(peak1.norm_scaled)*(peak2.norm_scaled)
                #pair_list is a list of tuples saying which two peaks are matched and the product of their intensities
                tuple = (peak1, peak2, product)
                peak_pairs.append(tuple)
    
    #sort matching peaks by descending intensity product
    peak_pairs.sort(reverse=True,key=lambda tuple: tuple[2])
    
    #keep list of which peaks have been used from each spectrum
    used_peaks_one=[]
    used_peaks_two=[]

    #reset total
    total = 0

    #loop through list of peak pairs in decreasing order of intensity
    for pair in peak_pairs:
        #split tuple
        peak_one=pair[0]
        peak_two=pair[1]
        product=pair[2]

        if (peak_one not in used_peaks_one) and (peak_two not in used_peaks_two):
            #sum intensity product if both peaks have not been used already
            total+=product
            
            #add to the list of used peaks
            used_peaks_one.append(peak_one)
            used_peaks_two.append(peak_two)
       
    return total

def main(file_path):
    """Give mgf file path as argument.
    Returns a list of spectrum pairs, showing the two spectrum IDs and the cosine similarity score of each pair"""


    #make list of spectrum objects from mgf file
    spectra_list=mgf_reader(file_path)

    spectra_matches=[]

    for spectrum_one in spectra_list:
        for spectrum_two in spectra_list:
            if spectrum_one.compare_precursor_mass(spectrum_two,1):#check if precursor masses are close enough to calculate cosine similarity
                
                #get spectrum IDs and calculate similarity
                spectrum_match=(spectrum_one.feature_id,spectrum_two.feature_id,cosine_score(spectrum_one,spectrum_two,0.3))
                spectra_matches.append(spectrum_match)
    
    return spectra_matches

if __name__ == "__main__":
    main(sys.argv[1])






