"""Methods to analyse MS1 intensities
"""

import csv
from scipy import stats

def sample_groups(groups_file):
    """Requires file with two columns 'samples', 'groups', determining which sample is in which of two groups.
    Returns a dictionary of the groups and the samples that belong in them
    """
    groups={}
    with open(groups_file) as csvfile:
        reader=csv.DictReader(csvfile)
        for row in reader:
            group=row['group']
            if group not in groups:
                groups[group]=[]
            groups[group].append(row['sample'])
    
    return groups


def samples_ttest(file_path,sample_groups_file,spectra_list):
    """Takes a csv file path of MS1 peak areas, a csv file path of which group each sample belongs in, and a list of Spectrum objects.
    Performs an independent t-test for each spectrum and adds the results to the Spectrum object.
    """
    groups=sample_groups(sample_groups_file)
    key1,key2 = groups.keys()

    group1=groups[key1]
    group2=groups[key2]

    with open(file_path) as csvfile:
        reader=csv.DictReader(csvfile)
        for row in reader:

            #peak areas of samples in first group
            list1=[]
            for sample in group1:
                result=[float(value) for key,value in row.items() if str(sample) in key and 'area' in key]
                if len(result)==1:
                    list1.append(result[0])
                elif len(result)==0:
                    print(f"No MS1 peak area found for sample {sample}")

            #peak areas of samples in second group
            list2=[]
            for sample in group2:
                result=[float(value) for key,value in row.items() if str(sample) in key and 'area' in key]
                if len(result)==1:
                    list2.append(result[0])
                elif len(result)==0:
                    print(f"No MS1 peak area found for sample {sample}")
            
            #t-test to compare the two groups
            statistic,p_value = stats.ttest_ind(list1,list2)
            #use p_value threshold of 0.05 for significance
            if p_value<=0.05:
                #determine which group has greater intensity value
                if statistic>0:
                    greater=key1
                else:
                    greater=key2
            else:
                greater="No significance"

            spectra_list.sort(key=lambda S: int(S.parameters['SCANS']))

            row_ID=[value for key,value in row.items() if 'ID' in key]
            row_ID=row_ID[0]

            #find spectrum that the row corresponds to and add t-test info
            added=False
            for i in range(len(spectra_list)):
                for S in spectra_list[i:]:
                    if S.parameters['SCANS'] == row_ID:
                        S.parameters['greater_group']=greater
                        S.parameters['tstatistic']=statistic
                        S.parameters['p_value']=p_value
                        added=True
                        break

                if added:
                    break




            


# samples_ttest("C:\\Users\\Kiah\\Documents\\Project\dummy_MS1.csv","C:\\Users\\Kiah\\Documents\\Project\\groups.csv")
# sample_groups("C:\\Users\\Kiah\\Documents\\Project\\groups.csv")
