"""Class to hold an MS2 peak
"""

import math

class Peak:
    """Class to hold an MS2 peak.

    Parameters:
    mass
    intensity
    sqrt_intensity -- calculated as soon as the peak is initialised
    scaled_intensity (after using Spectrum.euclidean_scale method) -- intensities in a spectrum are scaled to Euclidean norm 1
    """

    def __init__(self,mass,intensity):
        self.mass = mass
        self.intensity = intensity
        self.sqrt_intensity = math.sqrt(self.intensity)
    
    def __repr__(self):
        return f"{self.__class__.__name__}({self.mass},{self.intensity})"