# import math
# import numpy
# import sys
# import networkx as nx
# from functools import total_ordering
# import bisect
import network
import similarity
import read_mgf as mgf

def main(file_path):

    #make list of spectrum objects from mgf file
    print("reading file")
    spectra_list=mgf.read_mgf(file_path)[200:400]
   
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
    

##to use from command line by giving path to mgf file as argument
#if __name__ == "__main__":
#    main(sys.argv[1])


import time
start = time.time()

pairs=main(".\data\MS2_peaks.mgf")

elapsed = time.time()-start
#print start and end time of program running
print(time.strftime("\nelapsed time:\t%H:%M:%S", time.gmtime(elapsed)))




