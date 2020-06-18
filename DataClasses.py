import math
class Peak:
    def __init__(self,mass,intensity):
        self.mass = mass
        self.intensity = intensity
        self.sqrt_intensity = math.sqrt(self.intensity)
    
    def __repr__(self):
        return f"{self.__class__.__name__}({self.mass},{self.intensity})"

    def compare_peaks(self, peak_2, tolerance):
        if (self.mass >= peak_2.mass-tolerance and self.mass <= peak_2.mass+tolerance):
            return True
        else:
            return False

class Spectrum:
    def __init__(self):
        self.peaks = [] #list of peaks in spectrum
        
    def __repr__(self):
        return f"{self.__class__.__name__}({self.peaks})"

    def add_peak(self,mass, intensity):
        self.peaks.append(Peak(mass,intensity))
    
    def set_peak_id(self,peak_id):
        self.peak_id = peak_id

    def compare_precursor_mass(self,spectrum_2,threshold):
        print(f"{self.pep_mass}\t{spectrum_2.pep_mass}")
        if ((self.pep_mass-spectrum_2.pep_mass)>(-threshold) and (self.pep_mass-spectrum_2.pep_mass)<threshold):
            return True
        else:
            return False