import argparse
from mol_networking import similarity
import time
start = time.time()
'''
Usage:
    python molecular_networking.py .\data\MS2_peaks -o example -l .\massbank_library\MASSBANK
'''

parser = argparse.ArgumentParser(
    description='Test input arguments.',
    usage='%(prog)s <input-mgf> [options]',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)


parser.add_argument(
    'input',
    help='Input MGF file path and name, excluding \'.mgf\''
    )

parser.add_argument(
    '-o',
    '--output',
    help='Output GraphML file path and name, excluding \'.graphml\'',
    default='output'
    )

parser.add_argument(
    '-l',
    '--library',
    help='Library mgf file path and name, excluding \'.mgf\''
    )

parser.add_argument(
    '-lp',
    '--library-peaks',
    help='Minimum number of peaks required to match a library spectrum',
    type=int,
    default=3
    )

parser.add_argument(
    '-ls',
    '--library-score',
    help='Minimum similarity score required to match a library spectrum',
    type=float,
    default=0.7
    )

parser.add_argument(
    '-lpt',
    '--lib-precursor-tolerance',
    help='Precursor mass tolerance allowed when matching to library spectra',
    type=float,
    default=1.0
)

parser.add_argument(
    '--greedy',
    help='Use the fast greedy method of similarity scoring. Default will use maximum weighted method',
    action='store_true'
)

parser.add_argument(
    '-m',
    '--modified',
    help='Use modified cosine similarity scoring',
    default=True,
    type=bool
)

parser.add_argument(
    '-ma',
    '--modification-mass',
    help='Maximum modification mass allowed when calculating the modification value between two spectra',
    type=float,
    default=400.
)

parser.add_argument(
    '-ft',
    '--fragment-tolerance',
    help='Fragment tolerance when comparing two fragment peaks',
    type=float,
    default=0.1
)

parser.add_argument(
    '-p',
    '--peaks',
    help='Minimum number of peaks required to match a spectrum',
    type=int,
    default=6
)

parser.add_argument(
    '-s',
    '--score',
    help='Minumum cosine score for a spectral match to be kept',
    type=float,
    default=0.7
)

parser.add_argument(
    '-n',
    '--n-neighbours',
    help='Maximum number of neighbours a spectrum can have in the network',
    type=int,
    default=10
)

parser.add_argument(
    '-f',
    '--family-size',
    help='Maximum molecular family size allowed in the network',
    type=int,
    default=100
)

parser.add_argument(
    '--matchms',
    help='use the MatchMS for reading mgf file and calculating similarities',
    action='store_true'
)


args = parser.parse_args()
# print(args.accumulate(args.integers))
print(args)

if (args.matchms):
    print('usematchms')
    # import matchms
    from matchms.importing import load_from_mgf
    from matchms.filtering import default_filters
    from matchms.filtering import normalize_intensities
    from matchms import calculate_scores
    from matchms.similarity import ModifiedCosine
    from mol_networking.use_matchms import convert_matches as convert
    
    input_mgf=f'{args.input}.mgf'
    file=load_from_mgf(input_mgf)

    # Apply filters to clean and enhance each spectrum
    spectrums = []
    for spectrum in file:
        spectrum = default_filters(spectrum)
        # Scale peak intensities to maximum of 1
        spectrum = normalize_intensities(spectrum)
        spectrums.append(spectrum)

    if (args.library):
        from mol_networking.use_matchms import library_match as lib
        library_mgf=f'{args.library}.mgf'
        lib.library_match(spectrums,library_mgf,precursor_tol=args.lib_precursor_tolerance,cosine=args.library_score,n_peaks=args.library_peaks)
        

    scores = calculate_scores(references=spectrums,
                          queries=spectrums,
                          similarity_function=ModifiedCosine(tolerance=args.fragment_tolerance))


    spectra_matches=convert.convert_scores(scores)
    spectra_list=[]
    for s in spectrums:
        new = convert.convert_spectrum(s)
        print(new.parameters)
        spectra_list.append(new)
    print(spectra_list[4].feature_id)

    

else:
    from mol_networking import read_mgf as mgf

    input_mgf=f'{args.input}.mgf'
    print(f"reading file {input_mgf}")
    
    spectra_list=mgf.read_mgf(input_mgf)

    if (args.library):
        library_file=f'{args.library}.mgf'
        print("comparing to library")
        similarity.library_match(spectra_list,library_file,precursor_tol=args.lib_precursor_tolerance,cosine=args.library_score,n_peaks=args.library_peaks)

    #calculate modified cosines, comparing each spectrum to every other spectrum
    print("calculating cosine scores")
    spectra_matches=similarity.compare_all(spectra_list,fragment_tolerance=args.fragment_tolerance,modified=args.modified,precursor_tolerance=args.modification_mass,greedy=args.greedy)

    print(len(spectra_list))

#filter matches
print("filtering spectrum matches")
spectra_matches=similarity.filter_pairs(
    spectra_matches,
    cosine_threshold=args.score,
    peak_threshold=args.peaks)


from mol_networking import network    

print("making graph")
graph=network.make_network(spectra_list,spectra_matches)
        
print("filtering neighbours")
graph=network.filter_neighbors(graph,args.n_neighbours)

print("filtering family size")
graph=network.filter_family(graph,args.family_size)

output_file=f'{args.output}.graphml'
network.write_graphml(graph, output_file)
print(f"written to {output_file}")


elapsed = time.time()-start
#print time of program running
print(time.strftime("\nelapsed time:\t%H:%M:%S", time.gmtime(elapsed)))
print("\n")
