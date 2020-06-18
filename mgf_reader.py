from DataClasses import *

def mgf_reader(file_path):
    file = open(file_path)
    lines = file.read().split('\n')
    spectra_list=[]
    for line in lines[:100]:
        if 'BEGIN IONS' in line:
            #initiate new ion list to contain lines of the next ion
            ions = []
            #initiate new spectrum object
            spectrum = Spectrum()
            #print("begin ions")
            continue
   
        if 'END IONS' in line:
            #at end of current ion, parse the ion data
            print(f"full spectra: {ions}")

            for i in ions:
                if "=" in i:
                    parameter,value = i.split("=")
                    if "feature_id".casefold() in parameter.casefold():
                        spectrum.set_peak_id(value)
                    elif "pepmass".casefold() in parameter.casefold():
                        spectrum.pep_mass = float(value)
                    elif "scans".casefold() in parameter.casefold():
                        spectrum.scans = float(value)
                    elif "RTinseconds".casefold() in parameter.casefold():
                        spectrum.retention_time = float(value)
                    elif "charge".casefold() in parameter.casefold():
                        spectrum.charge = value
                    elif "mslevel".casefold() in parameter.casefold():
                        spectrum.ms_level = int(value)
                else:
                    mass, intensity = i.split()
                    mass = float(mass)
                    intensity = float(intensity)
                    spectrum.add_peak(mass, intensity)
            spectra_list.append(spectrum)  
         
        #add line to existing spectra info
        ions.append(line)
    file.close()
    return spectra_list

file_path = ".\example_data\MS2_peaks.mgf"
spectra_list=mgf_reader(file_path)

#for i in range(len(spectra_list)):
#    print(spectra_list[0].compare_precursor_mass(spectra_list[i],100))

def cosine_score(spectrum_one, spectrum_two, tolerance):
    if spectrum_one==spectrum_two:
        print("same")
        return
    
    pairs_list=[]
    for peak1 in spectrum_one.peaks:
        for peak2 in spectrum_two.peaks:
            if peak1.compare_peaks(peak2,tolerance):
                product=(peak1.sqrt_intensity)*(peak2.sqrt_intensity)
                tuple = (peak1, peak2, product)
                pairs_list.append(tuple)
    
    print(pairs_list)
    

    


print(f"spectra: {spectra_list}")
for spectrum_one in spectra_list:
    
    for spectrum_two in spectra_list:
        if spectrum_one.compare_precursor_mass(spectrum_two,400):
            print(f"True")
            cosine_score(spectrum_one,spectrum_two,0.3)

        else:
            print(f"False")


