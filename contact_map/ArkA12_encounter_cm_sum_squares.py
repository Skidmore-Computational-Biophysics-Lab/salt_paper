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
maxB = 72
traj_combined_list = []

# DEBUGGING, CHECK IF PRINT WORKS
print('print is working')

traj_no_salt_list = [
    "/cluster_0_allsim.mdcrd",
    "/cluster_1_allsim.mdcrd",
    "/cluster_2_allsim.mdcrd",
    "/cluster_3_allsim.mdcrd",
    "/cluster_4_allsim.mdcrd"]

prmtop_list = [
    "/data/chem_shared/SH3/ArkA12-SH3_binding/ArkA12_1_gas.prmtop",
    "/data/chem_shared/SH3/ArkA12-SH3_binding/ArkA12_1_gas.prmtop",
    "/data/chem_shared/SH3/ArkA12-SH3_binding/ArkA12_1_gas.prmtop",
    "/data/chem_shared/SH3/ArkA12-SH3_binding/ArkA12_1_gas.prmtop",
    "/data/chem_shared/SH3/ArkA12-SH3_binding/ArkA12_1_gas.prmtop"]

traj_no_salt_new = []
for traj_path, top_path in zip(traj_no_salt_list, prmtop_list):
    full_traj_path = "/data/mcohen3/SH3/no_salt_stripped_mdcrds" + traj_path
    traj_no_salt_new.append(
        md.load_netcdf(full_traj_path, top=top_path, stride=100) # the stride should be some multiple of 10
    )

arka12_stat = pd.read_csv('/data/mcohen3/SH3/no_salt_analysis/status_files_use_this/ArkA12_binding_no_salt_status.csv')
arka12_stat = arka12_stat.iloc[99::100, :]

num_clusters = 5
num_sims = 10

no_salt_split_traj = []

# split trajectories into their independent simulations
for j in range(len(traj_no_salt_new)):
    n = int(traj_no_salt_new[j].n_frames/10)
    
    split_traj = [traj_no_salt_new[j][i:i + n] for i in range(0, len(traj_no_salt_new[j]), n)]
    no_salt_split_traj.append(split_traj)

frames_per_sim = len(arka12_stat) // (num_clusters * num_sims)

traj_no_salt_encounter_ind = [] # for independent simulations

for i in range(num_clusters):
    for sim in range(num_sims):
        start = (i * num_sims + sim) * frames_per_sim
        stop = start + frames_per_sim

        arka12_stat_cluster_sim = arka12_stat.iloc[start:stop]

        mask = (
            (arka12_stat_cluster_sim['status'] != 'unbound') &
            (arka12_stat_cluster_sim['status'] != 'Bound')
        ).to_numpy()

        traj_no_salt_encounter_ind.append(no_salt_split_traj[i][sim][mask])

traj_no_salt_combined = md.join(traj_no_salt_encounter_ind)

traj_no_salt_combined_contacts = ContactFrequency(traj_no_salt_combined)

traj_salt_list = [
    "/cluster_0/cluster_0_allsim.mdcrd",
    "/cluster_1/cluster_1_allsim.mdcrd",
    "/cluster_2/cluster_2_allsim.mdcrd",
    "/cluster_3/cluster_3_allsim.mdcrd",
    "/cluster_4/cluster_4_allsim.mdcrd",]

prmtop_list = [
    "/data/chem_shared/SH3/ArkA12-SH3_binding_800mM_NaCl/cluster_0/ArkA12_0_gas_ions.prmtop",
    "/data/chem_shared/SH3/ArkA12-SH3_binding_800mM_NaCl/cluster_1/ArkA12_1_gas_ions.prmtop",
    "/data/chem_shared/SH3/ArkA12-SH3_binding_800mM_NaCl/cluster_2/ArkA12_2_gas_ions.prmtop",
    "/data/chem_shared/SH3/ArkA12-SH3_binding_800mM_NaCl/cluster_3/ArkA12_3_gas_ions.prmtop",
    "/data/chem_shared/SH3/ArkA12-SH3_binding_800mM_NaCl/cluster_4/ArkA12_4_gas_ions.prmtop",
]

traj_salt_new = []
for traj_path, top_path in zip(traj_salt_list, prmtop_list):
    full_traj_path = "/data/chem_shared/SH3/ArkA12-SH3_binding_800mM_NaCl" + traj_path
    traj_salt_new.append(
        md.load_netcdf(full_traj_path, top=top_path, stride=100) # the stride should be some multiple of 10
    )

arka12_stat = pd.read_csv('/data/jcardoso/binding_figures_800ArkA12/ArkA12/pair_wise_distances/binding/salt/status_files/ArkA12_binding_salt_status.csv')
arka12_stat = arka12_stat.iloc[99::100, :]

salt_split_traj = []

# split trajectories into their independent simulations
for j in range(len(traj_salt_new)):
    n = int(traj_salt_new[j].n_frames/num_sims)
    
    split_traj = [traj_salt_new[j][i:i + n] for i in range(0, len(traj_salt_new[j]), n)]
    salt_split_traj.append(split_traj)

num_clusters = 5
num_sims = 10

frames_per_sim = len(arka12_stat) // (num_clusters * num_sims)

traj_salt_encounter_ind = [] # for independent simulations

for i in range(num_clusters):
    for sim in range(num_sims):
        start = (i * num_sims + sim) * frames_per_sim
        stop = start + frames_per_sim

        arka12_stat_cluster_sim = arka12_stat.iloc[start:stop]

        mask = (
            (arka12_stat_cluster_sim['status'] != 'unbound') &
            (arka12_stat_cluster_sim['status'] != 'Bound')
        ).to_numpy()

        traj_salt_encounter_ind.append(salt_split_traj[i][sim][mask])

for traj in range(len(traj_salt_encounter_ind)):
    traj_salt_encounter_ind[traj] = traj_salt_encounter_ind[traj].remove_solvent()

traj_combined_salt = md.join(traj_salt_encounter_ind)

traj_combined_contacts_salt = ContactFrequency(traj_combined_salt)

diff = AtomMismatchedContactDifference(traj_no_salt_combined_contacts, traj_combined_contacts_salt)
difference = diff.residue_contacts.df.iloc[minB-1:maxB-1,minA-1:maxA]

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
    contacts_df = ContactFrequency(traj).residue_contacts.df.iloc[minB-1:maxB-1, minA-1:maxA]
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