from mol_networking import network, similarity
from mol_networking import read_mgf as mgf
import matchms

def main(file_path):

    #make list of spectrum objects from mgf file
    print("reading file")
    spectra_list=mgf.read_mgf(file_path)[200:300]
   
    #match to library file to add names to spectra objects
    print("comparing to library")
    similarity.library_match(spectra_list,".\massbank_library\MASSBANK.mgf")
    
    #calculate modified cosines, comparing each spectrum to every other spectrum
    print("calculating cosine scores")
    spectra_matches=similarity.compare_all(spectra_list,modified=True,precursor_tolerance=400)
        
    #filter matches
    print("filtering spectrum matches")
    spectra_matches=similarity.filter_pairs(spectra_matches)

    print("making graph")
    graph=network.make_network(spectra_list,spectra_matches)

    print("filtering neighbours")
    graph=network.filter_neighbors(graph,10)

    print("filtering family size")
    graph=network.filter_family(graph,100)

    network.write_graphml(graph, "subset_network.graphml")
    

import time
start = time.time()

pairs=main(".\data\MS2_peaks.mgf")

elapsed = time.time()-start
#print time of program running
print(time.strftime("\nelapsed time:\t%H:%M:%S", time.gmtime(elapsed)))