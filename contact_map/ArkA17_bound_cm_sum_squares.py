import mdtraj as md
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib
from contact_map import ContactFrequency, ContactDifference, AtomMismatchedContactDifference
import scipy.stats as st
from scipy.stats import ttest_ind
from tqdm.notebook import tqdm
import random

minA = 1
maxA = 58
minB = 60
maxB = 76
traj_combined_list = []

# DEBUGGING, CHECK IF PRINT WORKS
print('print is working')

traj_list = ["abp1-arka_1_stripped.mdcrd", "abp1-arka_2_stripped.mdcrd", "abp1-arka_3_stripped.mdcrd", "abp1-arka_4_stripped.mdcrd", "abp1-arka_5_stripped.mdcrd", 
          "abp1-arka_6_stripped.mdcrd", "abp1-arka_7_stripped.mdcrd", "abp1-arka_8_stripped.mdcrd", "abp1-arka_9_stripped.mdcrd", "abp1-arka_10_stripped.mdcrd"]

traj_no_salt_new = []
for suffix in tqdm(traj_list):
    traj_no_salt_new.append(md.load_netcdf("/data/oonwuzul/abp1-arka17_bound/no_salt/analysis/mdcrd/prod_stripped/" + suffix, 
                          top="/data/aadkins/SH3/abp1-arka17/abp1-arka_gas.prmtop", stride=300))

num_clusters = 1
num_sims = 10

no_salt_split_traj = []

# split trajectories into their independent simulations
for j in range(len(traj_no_salt_new)):
    n = int(traj_no_salt_new[j].n_frames/10)
    
    split_traj = [traj_no_salt_new[j][i:i + n] for i in range(0, len(traj_no_salt_new[j]), n)]
    no_salt_split_traj.append(split_traj)

frames_per_sim = 800

traj_no_salt_encounter_ind = [] # for independent simulations

for i in range(num_clusters):
    for sim in range(num_sims):
        start = (i * num_sims + sim) * frames_per_sim
        stop = start + frames_per_sim
        traj_no_salt_encounter_ind.append(no_salt_split_traj[i][sim])

traj_no_salt_combined = md.join(traj_no_salt_encounter_ind)

traj_no_salt_combined_contacts = ContactFrequency(traj_no_salt_combined)

traj_list_salt = ["abp1-arka_highsalt_stripped_with_ions_1.mdcrd", "abp1-arka_highsalt_stripped_with_ions_2.mdcrd", "abp1-arka_highsalt_stripped_with_ions_3.mdcrd", "abp1-arka_highsalt_stripped_with_ions_4.mdcrd", "abp1-arka_highsalt_stripped_with_ions_5.mdcrd", 
          "abp1-arka_highsalt_stripped_with_ions_6.mdcrd", "abp1-arka_highsalt_stripped_with_ions_7.mdcrd", "abp1-arka_highsalt_stripped_with_ions_8.mdcrd", "abp1-arka_highsalt_stripped_with_ions_9.mdcrd", "abp1-arka_highsalt_stripped_with_ions_10.mdcrd"]

traj_salt_new = []
for suffix in traj_list_salt:
    traj_salt_new.append(md.load_netcdf("/data/oonwuzul/abp1-arka17_bound/salt/abp1-arka17_800mM_salt/mdcrd/prod_stripped/" + suffix, 
                          top="/data/oonwuzul/abp1-arka17_bound/salt/abp1-arka17_800mM_salt/prmtop/stripped_with_ions.prmtop", stride=300 ))

salt_split_traj = []

# split trajectories into their independent simulations
for j in range(len(traj_salt_new)):
    n = int(traj_salt_new[j].n_frames/num_sims)
    
    split_traj = [traj_salt_new[j][i:i + n] for i in range(0, len(traj_salt_new[j]), n)]
    salt_split_traj.append(split_traj)

num_clusters = 1
num_sims = 10

frames_per_sim = 800

traj_salt_encounter_ind = [] # for independent simulations

for i in range(num_clusters):
    for sim in range(num_sims):
        start = (i * num_sims + sim) * frames_per_sim
        stop = start + frames_per_sim

        traj_salt_encounter_ind.append(salt_split_traj[i][sim])

for traj in range(len(traj_salt_encounter_ind)):
    traj_salt_encounter_ind[traj] = traj_salt_encounter_ind[traj].remove_solvent()

traj_combined_salt = md.join(traj_salt_encounter_ind)

traj_combined_contacts_salt = ContactFrequency(traj_combined_salt)

diff = AtomMismatchedContactDifference(traj_no_salt_combined_contacts, traj_combined_contacts_salt)
difference = diff.residue_contacts.df.iloc[minB-1:maxB,minA-1:maxA]

# DEBUGGING, CHECK IF MADE IT TO RANDOMIZATION TEST
print("made it to randomization test")

# Randomization test

test_stat = np.nansum((difference**2).to_numpy())
print("test statistic = ", test_stat)

# Pool together simulations from each set into one list

traj_all_ind = traj_no_salt_encounter_ind + traj_salt_encounter_ind

for traj in range(len(traj_all_ind)):
    traj_all_ind[traj] = traj_all_ind[traj].remove_solvent()


# Randomly select half of these sims for the first category (eg. no salt) and the other half in the second category (eg. salt)

# DO THIS 10,000 TIMES

tot_sims = int(len(traj_all_ind)/2)
# --- STEP 1: Pre-calculate outside the loop ---
# Compute contact matrices for individual trajectories once
precomputed_contacts = []
for traj in traj_all_ind:
    # Get dataframe/matrix for individual trajectory
    contacts_df = ContactFrequency(traj).residue_contacts.df.iloc[minB-1:maxB, minA-1:maxA]
    precomputed_contacts.append(contacts_df.to_numpy())

# Convert to a single 3D NumPy array: shape -> (num_trajectories, rows, cols)
contacts_array = np.array(precomputed_contacts)

# Compute test statistic using precomputed data
no_salt_avg = np.nanmean(contacts_array[:tot_sims], axis=0)
salt_avg = np.nanmean(contacts_array[tot_sims : 2 * tot_sims], axis=0)
diff_test = no_salt_avg - salt_avg
new_test_stat = np.nansum(diff_test ** 2)
print("new test statistic = ", new_test_stat)

# --- STEP 2: Fast Permutation Loop ---
random_stat_list = []
sample_size = 10000  # Example realistic sample size

for sample in range(sample_size):
    # Fast shuffle of indices instead of full objects
    shuffled_indices = np.random.permutation(len(traj_all_ind))
    
    group1_idx = shuffled_indices[:tot_sims]
    group2_idx = shuffled_indices[tot_sims : 2 * tot_sims]
    
    # Average precomputed matrices directly using NumPy
    avg_1 = np.nanmean(contacts_array[group1_idx], axis=0)
    avg_2 = np.nanmean(contacts_array[group2_idx], axis=0)
    
    # Calculate difference squared sum
    diff = avg_1 - avg_2
    random_stat = np.nansum(diff ** 2)
    random_stat_list.append(random_stat)

# --- STEP 3: Vectorized P-value ---
pval = sum(1 for x in random_stat_list if x > test_stat)/sample_size
new_pval = sum(1 for x in random_stat_list if x > new_test_stat)/sample_size
print(f"p-value = {pval}")
print(f"new p-value = {new_pval}")
#print(f"random_stat_list = {random_stat_list}")