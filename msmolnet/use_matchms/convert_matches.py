"""
Methods used to convert matchms Spectrum and Scores into MSMolNet objects to do molecular networking
"""

from msmolnet.Peak import Peak
from msmolnet.Spectrum import Spectrum

def convert_scores(scores):
    """Converts a matchms Scores object into a dictionary of matches to be used in molecular networking
    Returns: dict
    """
    matches={}
    for score in scores:

        (reference, query, score, n_matching) = score
        if reference == query:
            continue
        
        
        reference=reference.metadata['scans']
        
        query=query.metadata['scans']


        if reference not in matches:
            matches[reference]={}
        matches[reference][query]={'cosine':score,'peaks':n_matching}
        if query not in matches:
            matches[query]={}
        matches[query][reference]={'cosine':score,'peaks':n_matching}
        
    
    return matches


def convert_spectrum(matchms_spectrum):
    """Converts a matchms Spectrum object into a MSMolNet Spectrum object
    Returns: MSMolNet.Spectrum
    """
    peaks=matchms_spectrum.peaks
    mz=peaks.mz
    intensities=peaks.intensities
    spectrum=Spectrum()

    for i in range(len(mz)):
        spectrum.add_peak(mz[i],intensities[i])

    spectrum.parameters=matchms_spectrum.metadata

    if hasattr(matchms_spectrum,'library_parameters'):
        spectrum.library_parameters=matchms_spectrum.library_parameters

    if "pepmass" in spectrum.parameters:
        if  isinstance(spectrum.parameters['pepmass'],tuple):
            spectrum.parameters['pepmass']=spectrum.parameters['pepmass'][0]

    spectrum.set_id()
    return spectrum
