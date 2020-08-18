import math

class Peak:
    """Class to hold an MS2 peak"""
    def __init__(self,mass,intensity):
        self.mass = mass
        self.intensity = intensity
        self.sqrt_intensity = math.sqrt(self.intensity)
    
    def __repr__(self):
        return f"{self.__class__.__name__}({self.mass},{self.intensity})"
        # return str(self.mass)