import argparse

'''
Usage:
    python molecular_networking.py .\data\MS2_peaks -o example -l .\massbank_library\MASSBANK
'''

parser = argparse.ArgumentParser(
    description='Test input arguments.',
    usage='%(prog)s <input-mgf> [options]')


parser.add_argument(
    'input',
    help='input mgf file path and name, excluding \'.mgf\''
    )

parser.add_argument(
    '-o',
    '--output',
    help='output graphml name, excluding \'.graphml\'',
    default='output'
    )

parser.add_argument(
    '-l',
    '--library',
    help='library mgf file, excluding \'.mgf\''
    )

parser.add_argument(
    '-lp',
    '--library-peaks',
    help='number of peaks required to match a library spectrum',
    type=int,
    default=3
    )

parser.add_argument(
    '-lt',
    '--library-threshold',
    help='threshold similarity score required to match a library spectrum',
    type=float,
    default=0.7
    )

parser.add_argument(
    '-lpt',
    '--lib-precursor-tolerance',
    help='tolerance of precursor mass allowed when matching to library spectra',
    type=float,
    default=1.0
)

parser.add_argument(
    '--greedy',
    help='Use the fast greedy method of similarity scoring. Default will use maximum weighted method',
    action='store_true'
)

parser.add_argument(
    '--unmodified',
    help='use unmodified cosine (default: use modified cosine when calculating similarities)',
    action='store_true'
)

parser.add_argument(
    '-mt',
    '--modification-tolerance',
    help='m/z tolerance when calculating the modification value between two spectra',
    type=float,
    default=400.
)

parser.add_argument(
    '-ft',
    '--fragment-tolerance',
    help='m/z tolerance when comparing two fragment peaks',
    type=float,
    default=0.1
)

parser.add_argument(
    '-p',
    '--peaks',
    help='number of matching peaks required to match a spectrum',
    type=int,
    default=6
)

parser.add_argument(
    '-t',
    '--cosine-threshold',
    help='threshold cosine score for a spectral match to be kept',
    type=float,
    default=0.7
)

parser.add_argument(
    '-n',
    '--n-neighbours',
    help='maximum number of neighbours a spectrum can have in the network',
    type=int,
    default=10
)

parser.add_argument(
    '-f',
    '--family-size',
    help='maximum molecular family size allowed in the network',
    type=int,
    default=100
)

parser.add_argument(
    '--matchms',
    help='use the matchms package for reading mgf file and calculating similarities',
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
        lib.library_match(spectrums,library_mgf,precursor_tol=args.lib_precursor_tolerance,cosine=args.library_threshold,n_peaks=args.library_peaks)
        

    scores = calculate_scores(references=spectrums,
                          queries=spectrums,
                          similarity_function=ModifiedCosine(tolerance=args.fragment_tolerance))


    spectra_matches=convert.convert_scores(scores)
    spectra_list=[]
    for s in spectrums:
        spectra_list.append(convert.convert_spectrum(s))

    

else:
    from mol_networking import read_mgf as mgf
    from mol_networking import similarity

    input_mgf=f'{args.input}.mgf'
    print(f"reading file {input_mgf}")
    
    spectra_list=mgf.read_mgf(input_mgf)

    if (args.library):
        library_file=f'{args.library}.mgf'
        print("comparing to library")
        similarity.library_match(spectra_list,library_file,precursor_tol=args.lib_precursor_tolerance,cosine=args.library_threshold,n_peaks=args.library_peaks)

    #calculate modified cosines, comparing each spectrum to every other spectrum
    print("calculating cosine scores")
    spectra_matches=similarity.compare_all(spectra_list[:10],fragment_tolerance=args.fragment_tolerance,modified=not(args.unmodified),precursor_tolerance=args.modification_tolerance,greedy=args.greedy)

#filter matches
print("filtering spectrum matches")
spectra_matches=similarity.filter_pairs(
    spectra_matches,
    cosine_threshold=args.cosine_threshold,
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



