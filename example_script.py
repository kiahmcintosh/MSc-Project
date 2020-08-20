"""Example of MSMolNet workflow"""

from msmolnet import read_mgf, similarity, compare_ms1, network

#make list of spectrum objects from mgf file
input_file = 'example_file.mgf'
print(f'reading file {input_file}')
spectra_list=read_mgf.read_mgf(input_file)
   
#match to library spectra
print("comparing to library")
similarity.library_match(spectra_list,'library.mgf')

#ms1 analysis
print("MS1 t-test")
compare_ms1.samples_ttest('MS1_peaks.csv', 'groups.csv', spectra_list)
   
#calculate modified cosines for query spectra
print("calculating modified cosine scores")
spectra_matches=similarity.compare_all(spectra_list, modified=True, precursor_tolerance=400)
        
#filter matches using default thresholds
print("filtering spectrum matches")
spectra_matches=similarity.filter_pairs(spectra_matches)

#construct network
print("making network")
graph=network.make_network(spectra_list, spectra_matches)

#filter n neighbours
print("filtering n neighbours")
graph=network.filter_neighbors(graph, 8)

#filter molecular family size
print("filtering molecular family size")
graph=network.filter_family(graph, 50)

#write GraphML output file
output_file = 'network.graphml'
network.write_graphml(graph, output_file)
print(f"network written to {output_file}")