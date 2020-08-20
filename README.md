# MSMolNet

MSMolNet is a package used to carry out molecular networking. It contains a command line script to easily be able to run the molecular networking workflow.
<p>&nbsp;</p>

## Installation
Use git to clone the MSMolNet repository.

```
git clone https://github.com/kiahmcintosh/MSc-Project.git

```
Make a conda environment with the required packages.
```
conda create --name MSMolNet-env --file --MSMolNet-conda.spec
```
<p>&nbsp;</p>

## Command line script usage
Usage:
```
python run_MSMolNet <input-file> [options]
```

|Option|Default|Description|
|:-------------|:--------|:-------|
|\<input-file>||(Required) Input MGF file name, excluding '.mgf' 
|-h/--help|  | Display help options|
|-o/--output \<str> |'output'|File name prefix given to the output GraphML file|
|-p/--peaks \<int>|6|Minimum number of matching peaks required to match two spectra|
|-s/--score \<float>|0.7|Minimum similarity score required for a spectral match to be kept|
|--greedy||Use the fast greedy method of selecting peak matches instead of maximum-weighted|
|-m/--modified \<bool>|True|Use modified or unmodified cosine when calculating similarity|
|-ma/--modification-mass \<float>|400.0|Maximum modification mass allowed when comparing two spectra|
|-ft/--fragment-tolerance \<float>|0.3|Mass tolerance when comparing fragment peaks|
|-n/--n-neighbours \<int>|10|Maximum number of neighbours a spectrum can have in the network|
|-f/--family-size \<int>|100|Maximum molecular family size allowed in the network|
|-l/--library \<str>|None|MGF file name to be used for library matching, excluding '.mgf'|
|-ls/--library-score \<float>|0.7|Minimum cosine score required for a spectrum to match a library spectrum|
|-lp/--library-peaks \<int>|3|Minimum number of matching peaks required for a spectrum to match a library spectrum|
|-lpt/--lib-precursor-tolerance \<float>|1.0|Precursor mass tolerance for comparing a spectrum to library spectra|
|--matchms||Use matchms for reading the mgf file and calculating modified cosine scores|
|--ms1 \<str> \<str>||Carry out independent t-tests on MS1 feature intensities. Requires a CSV file of peak area in each sample for each spectrum and a CSV file with columns 'sample' and 'group'|
<p>&nbsp;</p>

## Command line examples

Create a network with file spectra_file.mgf and produce the output file network.graphml:
```
python run_MSMolNet.py spectra_file -o network
```
Alter some of the default filters for number of peaks, minimum score, and fragment mass tolerance for spectrum matching:
```
python run_MSMolNet.py spectra_file -o network -p 8 -s 0.8 -ft 0.2
```

Create a network and match the spectra to library spectra:
```
python run_MSMolNet.py spectra_file -o network -l library_file
```
Carry out MS1 analysis:
```
python run_MSMolNet.py spectra_file -o network --ms1 ms1_intensities.csv groups.csv
```
<p>&nbsp;</p>

## Script of example workflow
```
from msmolnet import read_mgf, similarity, compare_ms1, network

#make list of spectrum objects from mgf file
input_file = 'example_file.mgf'
print(f'reading file {input_file}')
spectra_list=read_mgf.read_mgf(input_file)
   
#match to library spectra
print("comparing to library")
similarity.library_match(spectra_list,'library.mgf')

#ms1 analysis
print("MS1 t-test")
compare_ms1.samples_ttest('MS1_peaks.csv', 'groups.csv', spectra_list)
   
#calculate modified cosines for query spectra
print("calculating modified cosine scores")
spectra_matches=similarity.compare_all(spectra_list, modified=True, precursor_tolerance=400)
        
#filter matches using default thresholds
print("filtering spectrum matches")
spectra_matches=similarity.filter_pairs(spectra_matches)

#construct network
print("making network")
graph=network.make_network(spectra_list, spectra_matches)

#filter n neighbours
print("filtering n neighbours")
graph=network.filter_neighbors(graph, 8)

#filter molecular family size
print("filtering molecular family size")
graph=network.filter_family(graph, 50)

#write GraphML output file
output_file = 'network.graphml'
network.write_graphml(graph, output_file)
print(f"network written to {output_file}")
```