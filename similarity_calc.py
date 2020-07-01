import math
import numpy
import sys
import networkx as nx
from functools import total_ordering
import bisect

class Peak:
    """Class to hold an MS2 peak"""
    def __init__(self,mass,intensity):
        self.mass = mass
        self.intensity = intensity
        self.sqrt_intensity = math.sqrt(self.intensity)
    
    def __repr__(self):
        return f"{self.__class__.__name__}({self.mass},{self.intensity})"
        
@total_ordering
class Spectrum:
    def __init__(self):
        self.peaks = [] #list of peaks in spectrum
        self.feature_id=None
        
    def __repr__(self):
        # return f"{self.__class__.__name__}({self.peaks})"
        return self.feature_id
    
    def __lt__(self,spectrum):
        if self.pep_mass<spectrum.pep_mass:
            return True
        else:
            return False


    def add_peak(self,mass, intensity):
        self.peaks.append(Peak(mass,intensity))
    
    def set_peak_id(self,feature_id):
        self.feature_id = feature_id

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

    # for line in lines[:250]:  ##alternative to test subset of data
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
                    parameter,value = i.split("=",1)
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
                    elif "name".casefold() in parameter.casefold():
                        spectrum.name=value

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

def calc_similarities(spectra_list,fragment_tolerance=0.3, modified=False,precursor_tolerance=1.0):
    """takes a list of spectrum objects and calculates modified cosine for each spectrum matched with every other spectrum.
    returns a list of spectrum matches with cosine score and number of matching peaks"""

    #list for matched spectra to be stored
    spectra_matches=[]
    
    for i in range(0,len(spectra_list)):
        spectrum_one=spectra_list[i]
        for j in range(i+1,len(spectra_list)):
            spectrum_two=spectra_list[j]
                            
            #cosine score calculation
            score,peak_count=cosine_score(spectrum_one,spectrum_two,fragment_tolerance,modified,precursor_tolerance)
                
            #if spectra don't match/precursor mass too far away, don't add to list
            if(score==0 and peak_count==0):
                continue

            #make tuple of spectra IDs,cosine score and number of matching peaks
            spectrum_match=(spectrum_one,spectrum_two,score,peak_count)
            # print(type(spectrum_match[0]))
            # print(type(spectrum_match[1]))
            # print(spectrum_match)
                
            spectra_matches.append(spectrum_match)
    
    return spectra_matches

def cosine_score(spectrum_one, spectrum_two, fragment_tolerance=0.3, modified=False,precursor_tolerance=1.0):
    """takes two MS2 Spectrum objects, with the fragment tolerance and returns the cosine similarity score.
    By default will calculate modified cosine; use modification=False to return normal cosine.
    Default fragment mass tolerance = 0.3
    Default precursor mass tolerance = 400"""
    modification=spectrum_two.pep_mass-spectrum_one.pep_mass
    if abs(modification)>precursor_tolerance:
            return 0,0
    
    if modified==False:
        modification=0
    
    peak_pairs=[]

    #make list of pairs of matched peaks within given tolerance
    for peak1 in spectrum_one.peaks:
        for peak2 in spectrum_two.peaks:

            #check normal and modified peak masses
            if (abs(peak1.mass-peak2.mass)<=fragment_tolerance) or (abs(peak1.mass+modification-peak2.mass) <= fragment_tolerance):
                product=(peak1.norm_scaled)*(peak2.norm_scaled)
                tuple = (peak1, peak2, product)
                #add to list of peak pairs
                peak_pairs.append(tuple)
            

    #sort matching peaks by descending intensity product
    peak_pairs.sort(reverse=True,key=lambda tuple: tuple[2])
    
    #keep list of which peaks have been used from each spectrum
    used_peaks_one=[]
    used_peaks_two=[]

    #reset total
    total = 0
    peak_count=0

    #loop through list of peak pairs in decreasing order of intensity
    for pair in peak_pairs:
        #split tuple
        peak_one=pair[0]
        peak_two=pair[1]
        product=pair[2]

        #keep track of number of peaks matched
        

        if (peak_one not in used_peaks_one) and (peak_two not in used_peaks_two):
            #sum intensity product if both peaks have not been used already
            total+=product
            #count number of matched peaks
            peak_count +=1
            
            #add to the list of used peaks
            used_peaks_one.append(peak_one)
            used_peaks_two.append(peak_two)

    return total, peak_count

def filter_pairs(pairs):
    """still to be completed
    Removes any matches with cosine <0.7 and <6 matches peaks"""
    filtered_pairs=[]
    for pair in pairs:
        #only keep spectra matches with cosine >0.7 and at least 6 matching peaks
        if pair[2]>=0.7 and pair[3]>=6:
            filtered_pairs.append(pair)
            
    # print("filtered:")
    # print(temp_pairs)
    # print(f"number of filtered: {len(filtered_pairs)}")


    #sort by cosine - will need to remove lowest scored to reduce family size
    # filtered_pairs.sort(key=lambda spectra_pair: spectra_pair[2])


    """rest of function to be added"""

    return filtered_pairs

def library_match(spectra_list,lib_mgf):
    """Reads a given library mgf file and matches the given spectra to the library spectra using normal cosine.
    Each test spectra is given the name of the library spectra match with the highest cosine score."""
    library=mgf_reader(lib_mgf)
    library_sort=sorted(library)

    for test_spectra in spectra_list:
        pos=bisect.bisect(library_sort,test_spectra)
        matches=[]
        for lib in library_sort[pos-2:pos+2]:
            score,peaks=cosine_score(test_spectra,lib,modified=False)
            if score>=0.7 and peaks>=3:
                matches.append((score,peaks,lib.name))
        
        if len(matches)>0:
            #sort possible library matches by cosine score
            matches.sort(reverse=True,key=lambda tuple: tuple[0])
            #use name of match with highest cosine score
            test_spectra.name=matches[0][2]
            # print(f"{test_spectra.feature_id}\t{test_spectra.name}")
            
def make_graph(nodes_list,edges_list):
    graph=nx.Graph()
    for node in nodes_list:
        graph.add_node(node)

    for edge in edges_list:
        graph.add_edge(edge[0],edge[1],weight=edge[2])
    
    print(nx.nodes(graph))

    nx.write_graphml(graph, "test.graphml")
    


def main(file_path):
    """Give mgf file path as argument.
    Returns a list of spectrum pairs, showing the two spectrum IDs, the cosine similarity score of each pair and number of matched peaks"""

    #make list of spectrum objects from mgf file
    print("reading file")
    spectra_list=mgf_reader(file_path)
   
    #match to library file to add names to spectra objects
    print("comparing to library")
    library_match(spectra_list,".\massbank_library\MASSBANK.mgf")
    
    #calculate modified cosines, comparing each spectrum to every other spectrum
    print("calculating cosine scores")
    spectra_matches=calc_similarities(spectra_list,modified=True,precursor_tolerance=400)

    #filter matches
    print("filtering spectrum matches")
    spectra_matches=filter_pairs(spectra_matches)

    print("making graph")
    make_graph(spectra_list,spectra_matches)
    
    # return spectra_matches

##to use from command line by giving path to mgf file as argument
#if __name__ == "__main__":
#    main(sys.argv[1])


import time
start = time.ctime()

pairs=main(".\example_data\MS2_peaks.mgf")

end = time.ctime()
#print start and end time of program running
print(f"start: {start}\nend: {end}")








