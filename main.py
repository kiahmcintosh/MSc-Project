from mol_networking import network, similarity
from mol_networking import read_mgf as mgf
from matchms.importing import load_from_mgf
from matchms.filtering import default_filters
from matchms.filtering import normalize_intensities
from matchms import calculate_scores
from matchms.similarity import CosineGreedy
from mol_networking.use_matchms import convert_matches as convert
import networkx as nx
from mol_networking import read_mgf as mgf
import matchms


def test_matchms():
    # Read spectrums from a MGF formatted file, for other formats see https://matchms.readthedocs.io/en/latest/api/matchms.importing.html
    # file = load_from_mgf(".\matchms\\tests\\pesticides.mgf")
    file = load_from_mgf(".\matchms\\tests\pesticides.mgf")


    # Apply filters to clean and enhance each spectrum
    spectrums = []
    for spectrum in file:
        # Apply default filter to standardize ion mode, correct charge and more.
        # Default filter is fully explained at https://matchms.readthedocs.io/en/latest/api/matchms.filtering.html .
        # spectrum = default_filters(spectrum)
        # Scale peak intensities to maximum of 1
        spectrum2 = normalize_intensities(spectrum)
        spectrums.append(spectrum2)

    print(spectrums[-1].peaks.intensities)

    # Calculate Cosine similarity scores between all spectrums
    # For other similarity score methods see https://matchms.readthedocs.io/en/latest/api/matchms.similarity.html .
    scores = calculate_scores(references=spectrums,
                          queries=spectrums,
                          similarity_function=CosineGreedy())

    # Print the calculated scores for each spectrum pair
    # for score in scores:

    #     (reference, query, score, n_matching) = score

    #     # Ignore scores between same spectrum and
    #     # pairs which have less than 20 peaks in common
    #     if reference != query and n_matching >= 20:
    #         print(f"Reference scan id: {reference.metadata['scans']}")
    #         print(f"Query scan id: {query.metadata['scans']}")
    #         print(f"Score: {score:.4f}")
    #         print(f"Number of matching peaks: {n_matching}")
    #         print("----------------------------")

   
   
    matches=convert.convert_scores(scores)
    nodes=[]
    for s in spectrums:
        nodes.append(convert.convert_spectrum(s))

    print("filtering spectrum matches")
    matches=similarity.filter_pairs(matches)
    print("making graph")
    graph=network.make_network(nodes,matches)
    # graph=nx.Graph(matches)
    # graph.add_nodes_from(nodes)

    # graph=nx.Graph(matches)
    # attributes={}
    # for N in nodes:
    #     graph.add_node(N)
    #     attributes[N]={}
    #     if hasattr(N,'parameters'):
    #         for P in N.parameters:
    #             print(f"{P}:\t{N.parameters[P]}")
    #             # attributes[N][P]=N.parameters[P]
    #     break
    # nx.set_node_attributes(graph,attributes)
    
    print("filtering neighbours")
    graph=network.filter_neighbors(graph,3)

    print("filtering family size")
    graph=network.filter_family(graph,5)

    output="pesticides_matchms_network.graphml"
    network.write_graphml(graph, output)
    print(f"written to {output}")


def main(file_path,output_file):

    #make list of spectrum objects from mgf file
    print("reading file")
    spectra_list=mgf.read_mgf(file_path)
    spectra_list=spectra_list[:30]
    print(len(spectra_list))
   
    #match to library file to add names to spectra objects
    # print("comparing to library")
    # similarity.library_match(spectra_list,".\massbank_library\MASSBANK.mgf")

    print("MS1 t-test")
    from mol_networking import compare_ms1
    compare_ms1.samples_ttest("C:\\Users\\Kiah\\Documents\\Project\dummy_MS1.csv","C:\\Users\\Kiah\\Documents\\Project\\groups.csv",spectra_list)
    
    #calculate modified cosines, comparing each spectrum to every other spectrum
    print("calculating cosine scores")
    spectra_matches=similarity.compare_all(spectra_list,modified=True,precursor_tolerance=400)
        
    # #filter matches
    print("filtering spectrum matches")
    spectra_matches=similarity.filter_pairs(spectra_matches)

    print("making graph")
    graph=network.make_network(spectra_list,spectra_matches)

    print("filtering neighbours")
    graph=network.filter_neighbors(graph,3)

    print("filtering family size")
    graph=network.filter_family(graph,5)

    network.write_graphml(graph, output_file)
    print(f"written to {output_file}")
    
def ms1(file_path,csv_file,groups_file):
    from mol_networking import compare_ms1
    print("reading file")
    spectra_list=mgf.read_mgf(file_path)

    compare_ms1.samples_ttest(csv_file,groups_file,spectra_list[:20])








import time
start = time.time()

# main(".\data\MS2_peaks.mgf","ms1_test.graphml")
# main(".\matchms\\tests\\pesticides.mgf","pesticides_my_network.graphml")
test_matchms()

# ms1(".\data\MS2_peaks.mgf","C:\\Users\\Kiah\\Documents\\Project\dummy_MS1.csv","C:\\Users\\Kiah\\Documents\\Project\\groups.csv")

elapsed = time.time()-start
#print time of program running
print(time.strftime("\nelapsed time:\t%H:%M:%S", time.gmtime(elapsed)))