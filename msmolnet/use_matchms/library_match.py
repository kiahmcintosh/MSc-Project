from matchms.similarity import CosineHungarian
from matchms.importing import load_from_mgf
from matchms.filtering import default_filters
from matchms.filtering import normalize_intensities
from matchms import calculate_scores

def library_match(spectra_list,lib_mgf,precursor_tol=1.0,cosine=0.7,n_peaks=3):
    """Reads a given library mgf file and matches the given spectra to the library spectra using normal cosine.
    Each test spectra is given the name of the library spectra match with the highest cosine score."""

    
    library=load_from_mgf(lib_mgf)

    # Apply filters to clean and enhance each spectrum
    library_spectra = []
    for spectrum in library:
        # spectrum = default_filters(spectrum)
        # Scale peak intensities to maximum of 1
        spectrum = normalize_intensities(spectrum)
        library_spectra.append(spectrum)


    scores = calculate_scores(references=library_spectra,
                          queries=spectra_list,
                          similarity_function=CosineHungarian())

    scores_list=[]
    for score in scores:
        print(score)
        scores_list.append(score)
    
    scores_list.sort(reverse=True,key=lambda tuple:tuple[2])


        
        



        # if reference != query and n_matching >= 20:

    # for test_spectra in spectra_list:
    #     pos=bisect.bisect(library_sort,test_spectra)
    #     matches=[]
    #     for lib in library_sort[pos-2:pos+2]:
    #         score,peaks=cosine_score_max(test_spectra,lib,modified=False,precursor_tolerance=precursor_tol)
    #         if score>=cosine and peaks>=n_peaks:
    #             matches.append((score,peaks,lib))
        
    #     if len(matches)>0:
    #         #sort possible library matches by cosine score
    #         matches.sort(reverse=True,key=lambda tuple: tuple[0])
    #         #use parameters of spectrum match with highest cosine score
    #         test_spectra.library_parameters=matches[0][2].parameters