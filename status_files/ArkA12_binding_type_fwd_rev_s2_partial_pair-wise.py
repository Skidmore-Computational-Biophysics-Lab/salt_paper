#This script categorizes an ArkA12-AbpSH3 simulation into binding categories and outputs a file called ArkA12_status.csv

#file_path - the path to the directory with all the .dat files that are made when using the contact_map.sh script
#pair_wise_file - the file with the pair-wise average distance by snapshot
#binding_distance_EC - the distance cut off for defining that ArkA is bound in the encounter complex
#binding_distance_FB - the distance cut off for defining that ArkA is fully bound
#the pdb you are using must be numbered with the SH3 domain first followed by ArkA
#this uses the first and last resiude of segment 1 to define forward or reverse binding



def type_binding(file_path, pair_wise_file, binding_distance_EC, binding_distance_FB):
 total = pd.read_csv(pair_wise_file, engine='python')
 
 # contacts associated with fwd interaction
 f1 = pd.read_table(file_path+"7-60.dat", sep='\s+')
 f1.columns = ['#Frame', '7-60']
 f2 = pd.read_table(file_path+"8-60.dat", sep='\s+')
 f2.columns = ['#Frame', '8-60']
 f3 = pd.read_table(file_path+"53-60.dat", sep='\s+')
 f3.columns = ['#Frame', '53-60']
 f4 = pd.read_table(file_path+"35-66.dat", sep='\s+')
 f4.columns = ['#Frame', '35-66']
 f5 = pd.read_table(file_path+"32-66.dat", sep='\s+')
 f5.columns = ['#Frame', '32-66']
 f6 = pd.read_table(file_path+"34-66.dat", sep='\s+')
 f6.columns = ['#Frame', '34-66']
 f7 = pd.read_table(file_path+"9-66.dat", sep='\s+')
 f7.columns = ['#Frame', '9-66']
 f8 = pd.read_table(file_path+"13-66.dat", sep='\s+')
 f8.columns = ['#Frame', '13-66']
 f9 = pd.read_table(file_path+"15-66.dat", sep='\s+')
 f9.columns = ['#Frame', '15-66']
 f10 = pd.read_table(file_path+"16-66.dat", sep='\s+')
 f10.columns = ['#Frame', '16-66']
 f11 = pd.read_table(file_path+"48-66.dat", sep='\s+')
 f11.columns = ['#Frame', '48-66']
 f12 = pd.read_table(file_path+"50-66.dat", sep='\s+')
 f12.columns = ['#Frame', '50-66']
 f52 = pd.read_table(file_path+"52-66.dat", sep='\s+')
 f52.columns = ['#Frame', '52-66']
 f53 = pd.read_table(file_path+"53-66.dat", sep='\s+')
 f53.columns = ['#Frame', '53-66']

 # contacts associated with reverse interaction
 r1 = pd.read_table(file_path+"7-66.dat", sep='\s+')
 r1.columns = ['#Frame', '7-66']
 r2 = pd.read_table(file_path+"8-66.dat", sep='\s+')
 r2.columns = ['#Frame', '8-66']
 r3 = pd.read_table(file_path+"53-66.dat", sep='\s+')
 r3.columns = ['#Frame', '53-66']
 r4 = pd.read_table(file_path+"35-60.dat", sep='\s+')
 r4.columns = ['#Frame', '35-60']
 r5 = pd.read_table(file_path+"32-60.dat", sep='\s+')
 r5.columns = ['#Frame', '32-60']
 r6 = pd.read_table(file_path+"34-60.dat", sep='\s+')
 r6.columns = ['#Frame', '34-60']

 # contacts associated with the c-terminus interacting in the correct orientation
 ct1 = pd.read_table(file_path+"13-71.dat", sep='\s+')
 ct1.columns = ['#Frame', '13-71']
 ct2 = pd.read_table(file_path+"14-71.dat", sep='\s+')
 ct2.columns = ['#Frame', '14-71']
 ct3 = pd.read_table(file_path+"15-71.dat", sep='\s+')
 ct3.columns = ['#Frame', '15-71']
 ct4 = pd.read_table(file_path+"16-71.dat", sep='\s+')
 ct4.columns = ['#Frame', '16-71']
 ct5 = pd.read_table(file_path+"48-71.dat", sep='\s+')
 ct5.columns = ['#Frame', '48-71']
 ct6 = pd.read_table(file_path+"31-70.dat", sep='\s+')
 ct6.columns = ['#Frame', '31-70']
 ct7 = pd.read_table(file_path+"32-70.dat", sep='\s+')
 ct7.columns = ['#Frame', '32-70']
 ct8 = pd.read_table(file_path+"35-70.dat", sep='\s+')
 ct8.columns = ['#Frame', '35-70']
 ct9 = pd.read_table(file_path+"48-70.dat", sep='\s+')
 ct9.columns = ['#Frame', '48-70']

 print("Finished reading in files!")

 binding = pd.concat([f1, total['avg'], f2['8-60'], f3['53-60'], f4['35-66'], f5['32-66'], f6['34-66'], r1['7-66'], r2 ['8-66'], r3['53-66'], r4['35-60'], r5['32-60'], r6['34-60'], ct1['13-71'], ct2['14-71'], ct3['15-71'],ct4['16-71'], ct5['48-71'], ct6['31-70'], ct7['32-70'], ct8['35-70'], ct9['48-70'], f7['9-66'], f8['13-66'], f9['15-66'], f10['16-66'],f11['48-66'], f12['50-66'], f52['52-66'] ], axis =1)

 conditions = [
    # unbound must have avg pairwise distance above a cutoff
    (binding['avg'] > float(binding_distance_EC)), 
    # fully bound must have avg pairwise distance below a cutoff
    (binding['avg'] < float(binding_distance_FB)), 
    # forward encounter must have K(3) to [8,9,10] and K(-3) to [33,35,36]
    (((binding['7-60']<8) | (binding['8-60']<8) | (binding['53-60']<8)) &((binding['35-66']<8)| (binding['32-66']<8)| (binding['34-66']<8))),
    # reverse encounter must have K(-3) to [8,9,10] and K(3) to [33,35,36]
    (((binding['7-66']<8) | (binding['8-66']<8) | (binding['53-66']<8)) &((binding['35-60']<8)| (binding['32-60']<8)| (binding['34-60']<8))), 
    # segment 2 only encounter must have K(-3) to [33,35,36] and [K(-8) to [14,15,16,17,49],L(-7) to [32,33,36,49]]
    (( (binding['35-66']<8) | (binding['32-66']<8) | (binding['34-66']<8) ) & ( ( (binding['13-71']<8) | (binding['14-71']<8) | (binding['15-71']<8) | (binding['16-71']<8) | (binding['48-71']<8) ) | ( (binding['31-70']<8) | (binding['32-70']<8) | (binding['35-70']<8) | (binding['48-70']<8) ) ))]
 
    
 choices = ["unbound","Bound","fwd","Rev", 'seg2 only']

 binding['status'] = np.select(conditions, choices, default='partial')

 print(binding['status'].value_counts(normalize=True)*100)
 status_only = binding[['#Frame', 'status']].copy()
 status_only.to_csv('ArkA12_binding_salt_status.csv')

 encounter = binding[ (binding['status'] != "unbound") & (binding['status'] != "Bound") ]
 print(encounter['status'].value_counts(normalize=True)*100)

if __name__ == '__main__':
 import sys
 import numpy as np
 import pandas as pd
 type_binding(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
