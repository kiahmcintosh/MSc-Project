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
        # return f"{self.__class__.__name__}({self.mass},{self.intensity})"
        return str(self.mass)
        
@total_ordering
class Spectrum:
    def __init__(self):
        self.peaks = [] #list of peaks in spectrum
        self.feature_id=None
        self.name=False
        
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

    # for line in lines[:5000]:  ##alternative to test subset of data
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
            # score,peak_count=cosine_score(spectrum_one,spectrum_two,fragment_tolerance,modified,precursor_tolerance)
            score,peak_count=max_cosine(spectrum_one,spectrum_two,fragment_tolerance,modified,precursor_tolerance)
                
            #if spectra don't match/precursor mass too far away, don't add to list
            if(score==0 and peak_count==0):
                continue

            #make tuple of spectra IDs,cosine score and number of matching peaks
            spectrum_match=(spectrum_one,spectrum_two,score,peak_count)
            
            spectra_matches.append(spectrum_match)
    
    return spectra_matches

def cosine_score(spectrum_one, spectrum_two, fragment_tolerance=0.3, modified=False,precursor_tolerance=1.0):
    """takes two MS2 Spectrum objects and returns the cosine similarity score and number of matched peaks.
    By default will calculate normal cosine; use modification=True to return modified cosine.
    Default fragment mass tolerance = 0.3
    Default precursor mass tolerance = 1.0"""
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

        if (peak_one not in used_peaks_one) and (peak_two not in used_peaks_two):
            #sum intensity product if both peaks have not been used already
            total+=product
            #count number of matched peaks
            peak_count +=1
            
            #add to the list of used peaks
            used_peaks_one.append(peak_one)
            used_peaks_two.append(peak_two)
    
    return total, peak_count

def max_cosine(spectrum_one, spectrum_two, fragment_tolerance=0.3, modified=False,precursor_tolerance=1.0):
    """takes two MS2 Spectrum objects and returns the maximum cosine similarity score and number of matched peaks.
    By default will calculate normal cosine; use modification=True to return modified cosine.
    Default fragment mass tolerance = 0.3
    Default precursor mass tolerance = 1.0"""

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
    
    #make graph containing the peaks in the two spectra
    B=nx.Graph()

    #add edges to graph
    for p in peak_pairs:
        B.add_edge(p[0],p[1],weight=p[2])    
    
    #maximum weighted score
    matching=nx.algorithms.max_weight_matching(B)

    #number of peaks matched between the two spectra
    peak_count=len(matching)

    total=0
    #sum the intensity products for the matched peaks
    for m in matching:
        for p in peak_pairs:
            if (m[0]==p[0] and m[1]==p[1]) or (m[0]==p[1] and m[1]==p[0]):
                total+=p[2]

    return total,peak_count
      
def filter_pairs(pairs):
    """Removes any spectra matches with cosine <0.7 and <6 matched peaks"""
    filtered_pairs=[]
    for pair in pairs:
        #only keep spectra matches with cosine >0.7 and at least 6 matching peaks
        if pair[2]>=0.7 and pair[3]>=6:
            filtered_pairs.append(pair)
            
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
            
def make_graph(nodes_list,edges_list):
    graph=nx.Graph()
    for node in nodes_list:
        if node.name:
            graph.add_node(node,colour=1,ID=node.name)
        else:
            graph.add_node(node,colour=0)
        

    for edge in edges_list:
        graph.add_edge(edge[0],edge[1],weight=edge[2])

    return graph

def filter_neighbours(spectra_list,matches,n):
    """Pass a list of Spectrum objects, a list of spectral matching tuples and n, the maximum number of neighbours (matches) each spectrum
    can have. If a spectrum has more than n neighbours, the matches with the lowest cosine score are removed until there are only n remaining.
    Returns the filtered list of spectral matches"""
    for spectrum in spectra_list:
        #new list to store neighbours
        neighbours=[]

        for m in matches:
            if spectrum == m[0] or spectrum == m[1]:
                neighbours.append(m)
        if len(neighbours)>n: #number of neighbours higher than given threshold
            #sort descending
            neighbours.sort(reverse=True,key=lambda tuple: tuple[2])

            #remove spectral matches so matches up to n are kept
            for R in neighbours[n:]:
                matches.remove(R)
            
    return matches
            
def get_family(matches,spectrum,family_list,family_matches):
    if spectrum not in family_list:
        family_list.append(spectrum)
    
    neighbours=[]
    for m in matches:
        if spectrum == m[0]:
            neighbours.append(m[1])
            if m not in family_matches:
                family_matches.append(m)

        elif spectrum == m[1]:
            neighbours.append(m[0])
            if m not in family_matches:
                family_matches.append(m)

    for n in neighbours:
        if n not in family_list:
            family_list.append(n)
            family_list,family_matches=get_family(matches,n,family_list,family_matches)
    
    return family_list,family_matches

def filter_family(matches,spectra,n):
    used_spectra=[]
    for s in spectra:
        if s not in used_spectra:
            while True:

                family,f_pairs=get_family(matches,s,[],[])

                if len(family)<=n:
                    break

                f_pairs.sort(key=lambda tuple: tuple[2])
                matches.remove(f_pairs[0])
            
            for f in family:
                used_spectra.append(f)
    
    return matches



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

    print("filtering neighbours")
    spectra_matches=filter_neighbours(spectra_list,spectra_matches,10)

    print("filtering family size")
    spectra_matches=filter_family(spectra_matches,spectra_list,100)

    print("making graph")
    graph=make_graph(spectra_list,spectra_matches)

    nx.write_graphml(graph, "filtered_graph.graphml")
    return spectra_matches

##to use from command line by giving path to mgf file as argument
#if __name__ == "__main__":
#    main(sys.argv[1])


import time
start = time.ctime()

pairs=main(".\example_data\MS2_peaks.mgf")

end = time.ctime()
#print start and end time of program running
print(f"start: {start}\nend: {end}")


