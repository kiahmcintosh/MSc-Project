from .Spectrum import Spectrum

def read_mgf(file_path):
    """Takes the path to a .mgf file. Returns a list of Spectrum objects, each containing the attributes contained in the
    file and a list of Peak objects"""
    file = open(file_path)
    #split file into list of lines
    lines = file.read().split('\n')

    #initiate list to hold all the Spectrum objects
    spectra_list=[]

    for line in lines:
        if 'BEGIN IONS' in line:
            spectrum_lines = [] #initiate new spectra list to contain lines of the next spectra data
            spectrum = Spectrum() #initiate new spectrum object
            continue
   
        if 'END IONS' in line:
            #at end of current spectra, parse the data
            
            
            for i in spectrum_lines:
                
                #extract parameter values and assign to attributes
                if "=" in i:
                    parameter,value = i.split("=",1)
                    if "feature_id".casefold() in parameter.casefold():
                        spectrum.set_peak_id(value)
                    elif "pepmass".casefold() in parameter.casefold():
                        spectrum.pep_mass = float(value)
                    # elif "scans".casefold() in parameter.casefold():
                    #     spectrum.scans = float(value)
                    # elif "RTinseconds".casefold() in parameter.casefold():
                    #     spectrum.retention_time = float(value)
                    # elif "charge".casefold() in parameter.casefold():
                    #     spectrum.charge = value
                    # elif "mslevel".casefold() in parameter.casefold():
                    #     spectrum.ms_level = int(value)
                    elif "name".casefold() in parameter.casefold():
                        spectrum.name=value

                    spectrum.parameters[parameter]=value

                #when it reaches a line with no '=' then it is one of the peaks, which is added to the spectrum
                else:
                    mass, intensity = i.split()
                    mass = float(mass)
                    intensity = float(intensity)
                    spectrum.add_peak(mass, intensity)

            #intensities scaled by Euclidean norm
            spectrum.euclidean_scale()

            spectra_list.append(spectrum)  
         
        #add line to existing spectra info
        spectrum_lines.append(line)

    file.close()
    return spectra_list
