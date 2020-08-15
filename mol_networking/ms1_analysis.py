import csv
from scipy import stats

def sample_groups(groups_file):
    groups={}
    with open(groups_file) as csvfile:
        reader=csv.DictReader(csvfile)
        for row in reader:
            group=row['group']
            if group not in groups:
                groups[group]=[]
            groups[group].append(row['sample'])
    
    return groups


def sig_change(a,b):
    statistic,p_value = stats.ttest_ind(a,b)
    if p_value<=0.05:
        if result[0]>1:
            return 1
        else:
            return (-1)
    else:
        # print("groups are not significantly different")
        return 0


def samples_ttest(file_path,sample_groups_file,spectra_list):
    # print(sample_groups(sample_groups_file))
    groups=sample_groups(sample_groups_file)
    key1,key2 = groups.keys()

    group1=groups[key1]
    group2=groups[key2]

    with open(file_path) as csvfile:
        reader=csv.DictReader(csvfile)
        for row in reader:

            list1=[]
            for sample in group1:
                result=[float(value) for key,value in row.items() if str(sample) in key and 'area' in key]
                if len(result)==1:
                    list1.append(result[0])
                elif len(result)==0:
                    print(f"No MS1 peak area found for sample {sample}")

            list2=[]
            for sample in group2:
                result=[float(value) for key,value in row.items() if str(sample) in key and 'area' in key]
                if len(result)==1:
                    list2.append(result[0])
                elif len(result)==0:
                    print(f"No MS1 peak area found for sample {sample}")
            
            statistic,p_value = stats.ttest_ind(list1,list2)
            if p_value<=0.05:
                if statistic>0:
                    greater=key1
                else:
                    greater=key2
            else:
                greater="No significance"
            
            # result=sig_change(list1,list2)
            # if result>0:
            #     greater=key1
            # elif result<0:
            #     greater=key2
            # else:
            #     greater="No significance"

            spectra_list.sort(key=lambda S: int(S.parameters['SCANS']))
            # print(spectra_list)

            row_ID=[value for key,value in row.items() if 'ID' in key]
            row_ID=row_ID[0]

            added=False
            for i in range(len(spectra_list)):
                # print(f"i: {i}")
                for S in spectra_list[i:]:
                    if S.parameters['SCANS'] == row_ID:
                        # print(f"row: {row_ID}\tspectrum: {S.parameters['SCANS']}")
                        S.parameters['greater_group']=greater
                        S.parameters['tstatistic']=statistic
                        S.parameters['p_value']=p_value
                        # print(f"sample: {S}\tgreater: {greater}")
                        added=True
                        break

                if added:
                    break




            


# samples_ttest("C:\\Users\\Kiah\\Documents\\Project\dummy_MS1.csv","C:\\Users\\Kiah\\Documents\\Project\\groups.csv")
# sample_groups("C:\\Users\\Kiah\\Documents\\Project\\groups.csv")
