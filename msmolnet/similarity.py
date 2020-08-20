"""
Methods to carry out cosine similarity calculations and filtering, and library matching
"""

import networkx as nx
import bisect
from msmolnet import network
from msmolnet.Peak import Peak
from msmolnet.Spectrum import Spectrum
from msmolnet import read_mgf as mgf


def compare_all(spectra_list,fragment_tolerance=0.3, modified=False,precursor_tolerance=1.0,greedy=False):
    """Takes a list of spectrum objects and calculates modified cosine for each spectrum matched with every other spectrum.
    returns a nested dictionary of spectrum matches with cosine score and number of matching peaks"""

    #dictionary to store spectral matches
    matches = {}   
    
    for i in range(0,len(spectra_list)):
        
        # # optional print of list position to track progress
        # print(i)

        spectrum_one=spectra_list[i]
        for j in range(i+1,len(spectra_list)):
            spectrum_two=spectra_list[j]
                            
            #cosine score calculation
            if greedy:
                score,peak_count=cosine_score_greedy(spectrum_one,spectrum_two,fragment_tolerance,modified,precursor_tolerance)
            else:
                score,peak_count=cosine_score_max(spectrum_one,spectrum_two,fragment_tolerance,modified,precursor_tolerance)
    
            #if spectra don't match/precursor mass too far away, don't add to list
            if(score==0 and peak_count==0):
                continue

            #add spectral match to dictionary
            if spectrum_one not in matches:
                matches[spectrum_one]={}
            matches[spectrum_one][spectrum_two]={'cosine':score,'peaks':peak_count}

            if spectrum_two not in matches:
                matches[spectrum_two]={}
            matches[spectrum_two][spectrum_one]={'cosine':score,'peaks':peak_count}
            
    return matches

def cosine_score_greedy(spectrum_one, spectrum_two, fragment_tolerance=0.3, modified=False,precursor_tolerance=1.0):
    """Takes two MS2 Spectrum objects and returns the cosine similarity score and number of matched peaks.
    By default will calculate normal cosine; use modification=True to return modified cosine scores.
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
                product=(peak1.scaled_intensity)*(peak2.scaled_intensity)
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

def cosine_score_max(spectrum_one, spectrum_two, fragment_tolerance=0.3, modified=False,precursor_tolerance=1.0):
    """Takes two MS2 Spectrum objects and returns the maximum cosine similarity score and number of matched peaks.
    Uses maximum-weighted method to select peak matches.
       
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
                product=(peak1.scaled_intensity)*(peak2.scaled_intensity)
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
      
def filter_pairs(pairs,cosine_threshold=0.7,peak_threshold=6):
    """Takes a dictionary of spectra matches and removes any matches below given cosine and n peaks thresholds
    Default score threshold = 0.7
    Default peaks threshold = 6
    """
    for spectrum in pairs:
        
        #only keep spectral matches if they pass the thresholds
        filtered = {key:val for key, val in pairs[spectrum].items() if val['cosine']>=cosine_threshold and val['peaks']>=peak_threshold}
        pairs[spectrum]=filtered
    
   
    #new dictionary to hold filtered matches
    filtered_pairs={}
    #only keep spectra with matches
    for spectrum in pairs:
        if bool(pairs[spectrum]):
            filtered_pairs[spectrum]=pairs[spectrum]
            
    return filtered_pairs

def library_match(spectra_list,lib_mgf,precursor_tol=1.0,cosine=0.7,n_peaks=3):
    """Reads a given library mgf file and matches the given spectra to the library spectra using normal cosine.
    Each test spectra is given the name of the library spectra match with the highest cosine score."""
    library=mgf.read_mgf(lib_mgf)
    library_sort=sorted(library)

    for test_spectra in spectra_list:
        #find position where test spectrum would go in library list by precursor mass
        pos=bisect.bisect(library_sort,test_spectra)
        matches=[]
        #calculate scores for positions above and below insertion point
        for lib in library_sort[pos-2:pos+2]:
            score,peaks=cosine_score_max(test_spectra,lib,modified=False,precursor_tolerance=precursor_tol)
            if score>=cosine and peaks>=n_peaks:
                matches.append((score,peaks,lib))
        
        if len(matches)>0:
            #sort possible library matches by cosine score
            matches.sort(reverse=True,key=lambda tuple: tuple[0])
            #use parameters of spectrum match with highest cosine score
            test_spectra.library_parameters=matches[0][2].parameters
            
