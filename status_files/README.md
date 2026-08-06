# How to create status file
1. Run the contact_distance_all.ptraj script to create .dat files that contain the distances between specific CA atoms of the ArkA peptide and SH3 domain. For different protein-protein interactions, this script can be modified to target residue interactions that are deemed "most important" during binding. This script contains contact distances for every residue.
2. Run the binding_distances_csv.ipynb notebook to create a csv file that will be used in the final status file generation, using only the important contact ditance .dat files.
3. Submit the .scr script that runs the py script, or run the following line in a folder that contains all of the aforementioned filetypes:
   \# python [binding_type_fwd_rev_s2_partial_pair-wise.py] [CONTACT DISTANCE .dat FILES PATHWAY] [.csv FILE GENERATED FROM .ipynb] [UNBOUND STATE CUTOFF] [BOUND STATE CUTOFF]
   python ArkA12_binding_type_fwd_rev_s2_partial_pair-wise.py distance_files/ salt_binding_distances.csv 23 11.5
