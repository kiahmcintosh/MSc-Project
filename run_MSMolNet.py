'''
Script to run MSMolNet from command line


Example usage:
    python run_MsMolNet.py <input_file> -o <output_file> -l <library_file>


usage: run_MSMolNet.py <input-mgf> [options]

Test input arguments.

positional arguments:
  input                 Input MGF file path and name, excluding '.mgf'

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Output GraphML file path and name, excluding
                        '.graphml' (default: output)
  -p PEAKS, --peaks PEAKS
                        Minimum number of peaks required to match a spectrum
                        (default: 6)
  -s SCORE, --score SCORE
                        Minumum cosine score for a spectral match to be kept
                        (default: 0.7)
  --greedy              Use the fast greedy method of similarity scoring.
                        Default will use maximum weighted method (default:
                        False)
  -m MODIFIED, --modified MODIFIED
                        Use modified cosine similarity scoring (default: True)
  -ma MODIFICATION_MASS, --modification-mass MODIFICATION_MASS
                        Maximum modification mass allowed when calculating the
                        modification value between two spectra (default:
                        400.0)
  -ft FRAGMENT_TOLERANCE, --fragment-tolerance FRAGMENT_TOLERANCE
                        Fragment tolerance when comparing two fragment peaks
                        (default: 0.3)
  -n N_NEIGHBOURS, --n-neighbours N_NEIGHBOURS
                        Maximum number of neighbours a spectrum can have in
                        the network (default: 10)
  -f FAMILY_SIZE, --family-size FAMILY_SIZE
                        Maximum molecular family size allowed in the network
                        (default: 100)
  -l LIBRARY, --library LIBRARY
                        Library mgf file path and name, excluding '.mgf'
                        (default: None)
  -ls LIBRARY_SCORE, --library-score LIBRARY_SCORE
                        Minimum similarity score required to match a library
                        spectrum (default: 0.7)
  -lp LIBRARY_PEAKS, --library-peaks LIBRARY_PEAKS
                        Minimum number of peaks required to match a library
                        spectrum (default: 3)
  -lpt LIB_PRECURSOR_TOLERANCE, --lib-precursor-tolerance LIB_PRECURSOR_TOLERANCE
                        Precursor mass tolerance allowed when matching to
                        library spectra (default: 1.0)
  --matchms             use the MatchMS for reading mgf file and calculating
                        similarities (default: False)
  --ms1 MS1 MS1         Do t-test on MS1 data. Requires a .csv file of peak
                        area in each sample for each spectrum and a .csv file
                        with columns "sample" and "group" (default: None)

'''

import argparse
from msmolnet import similarity

# #used to measure run time
# import time
# start = time.time()


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
    default=0.3
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
    '-l',
    '--library',
    help='Library mgf file path and name, excluding \'.mgf\''
    )

parser.add_argument(
    '-ls',
    '--library-score',
    help='Minimum similarity score required to match a library spectrum',
    type=float,
    default=0.7
    )

parser.add_argument(
    '-lp',
    '--library-peaks',
    help='Minimum number of peaks required to match a library spectrum',
    type=int,
    default=3
    )

parser.add_argument(
    '-lpt',
    '--lib-precursor-tolerance',
    help='Precursor mass tolerance allowed when matching to library spectra',
    type=float,
    default=1.0
)


parser.add_argument(
    '--matchms',
    help='use the MatchMS for reading mgf file and calculating similarities',
    action='store_true'
)

parser.add_argument(
    '--ms1',
    help='''Do t-test on MS1 data. Requires a .csv file of peak area in each sample for each spectrum 
    and a .csv file with columns "sample" and "group"''',
    nargs=2
    
)


args = parser.parse_args()
print(args)

if (args.matchms):
    print('use matchms')
    import matchms
    from matchms.importing import load_from_mgf
    from matchms.filtering import default_filters
    from matchms.filtering import normalize_intensities
    from matchms.filtering import add_precursor_mz
    from matchms import calculate_scores
    from matchms.similarity import ModifiedCosine
    from msmolnet.use_matchms import convert_matches as convert
    
    input_mgf=f'{args.input}.mgf'
    print(f"reading file {input_mgf}")
    file=load_from_mgf(input_mgf)
    print(file)

    print("normalising intensities")
    
    # Apply filters to clean and enhance each spectrum
    spectrums = []
    for spectrum in file:
        spectrum = default_filters(spectrum)
        # Scale peak intensities to maximum of 1
        spectrum = normalize_intensities(spectrum)
        print(spectrum.get('precursor_mz'))
        spectrums.append(spectrum)
 

    scores = calculate_scores(references=spectrums,
                          queries=spectrums,
                          similarity_function=ModifiedCosine(tolerance=args.fragment_tolerance))


    spectra_matches=convert.convert_scores(scores)
    spectra_list=[]
    for s in spectrums:
        new = convert.convert_spectrum(s)
        spectra_list.append(new)


    

else:
    from msmolnet import read_mgf as mgf

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


if (args.ms1):
    from msmolnet import ms1_analysis
    ms1_analysis.samples_ttest(args.ms1[0],args.ms1[1],spectra_list)

#filter matches
print("filtering spectrum matches")
spectra_matches=similarity.filter_pairs(
    spectra_matches,
    cosine_threshold=args.score,
    peak_threshold=args.peaks)


from msmolnet import network    

print("making graph")
graph=network.make_network(spectra_list,spectra_matches)
        
print("filtering neighbours")
graph=network.filter_neighbors(graph,args.n_neighbours)

print("filtering family size")
graph=network.filter_family(graph,args.family_size)

output_file=f'{args.output}.graphml'
network.write_graphml(graph, output_file)
print(f"written to {output_file}")

# #print time of program running
# elapsed = time.time()-start
# print(time.strftime("\nelapsed time:\t%H:%M:%S", time.gmtime(elapsed)))
# print("\n")
