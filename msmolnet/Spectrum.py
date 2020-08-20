"""Class to hold an MS2 spectrum
"""

import numpy
from .Peak import Peak
from functools import total_ordering

@total_ordering
class Spectrum:
    """Class to hold an MS2 Spectrum

    Parameters:
    peaks -- a list of Peak objects
    feature_id -- set from 'scans' in metadata, used as node labels in the network
    parameters -- dictionary of metadata provided in the MGF file
    library_parameters (if similarity.library_match method has been used) --metadata from the matched library spectrum

    """
    def __init__(self):
        self.peaks = [] #list of peaks in spectrum
        self.feature_id=None
        # self.name=False
        self.parameters={}
        
    def __repr__(self):
        # return f"{self.__class__.__name__}({self.peaks})"
        return self.feature_id
    
    def __lt__(self,spectrum):
        if self.pep_mass<spectrum.pep_mass:
            return True
        else:
            return False


    def add_peak(self,mass, intensity):
        """Create a peak object and add to the list of spectrum peaks
        """
        self.peaks.append(Peak(mass,intensity))

    def euclidean_scale(self):
        """Scales the intensity of each Peak object using Euclidean norm"""
        
        #get list of intensities from spectrum
        intensities=[]
        for peak in self.peaks:
            intensities.append(peak.sqrt_intensity)
        
        #calculate norm
        norm=numpy.linalg.norm(intensities)

        #set each peak's scaled intensity
        for peak in self.peaks:
            peak.scaled_intensity=(peak.sqrt_intensity)/norm

    def set_id(self):
        """Get scans number from metadata and set as spectrum ID
        """
        if 'SCANS' in self.parameters:
            self.feature_id=self.parameters['SCANS']
        if 'scans' in self.parameters:
            self.feature_id=self.parameters['scans']

    def __str__(self):
        return str(self.feature_id)