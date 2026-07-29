import pandas as pd
import numpy as np
import math
from scipy import stats
import matplotlib.pyplot as plt


# function to calculate ion interactions with charged residues
# for binding simulaitons, only includes the frames in the encounter complex
# output:
#   - full_ion_contact_data_label.csv - fraction of ensemble in contact with an ion of opposite charge for each residue
#   - full_avg_ion_data_label.csv - average number of ions in contact with residue at any point in time
#   - full_avg_time_diff_ion_data_label.csv - average total time that a particular ion stays in contact with that residue during one simulation in ns (this will be longer for longer simulations, simulations with less total ions present)
#   - full_avg_num_diff_ion_data_label.csv - average number of distinct ions that are in contact with that residue in one simulation (this will depend on the number of ions in the box, simulation length, as well as the frequency of interactions)
#   - full_total_contacts_data_label.csv - average total number of residues in contact with cations or anions at any given time
#   - full_total_ions_data_label.csv - average total number of anions or cations in contact with a charged residue at any given time
#   - figures/800mM_NaCl_ion_contacts_label.pdf - plot of fraction of ensemble in contact with an ion of opposite charge for each residue
def ion_contact(filename, residue, min, max, nsim, statusfile):
    # filename is the path and beginning of the data file names
    # residue is the protein residue that we are looking for ion contacts with (must be a string)
    # min and max are the first and last residue numbers corresponding to the salt ions
    # nsim is the number of independent simulations run (code assumes all simulations are the same length)
    
    # read in status for binding simulations
    if statusfile != "none":
        stat = pd.read_csv(statusfile)

    files_list = [] # list of dataframes for each ion
    sim_averages = [] # list that will contain the fraction of ensemble in contact with an ion for each independent simulation
    sim_ions_averages = [] # list that will contain the average number of ions in contact
    # code to calc for each distinct ion
    diff_ions_averages = [] # list that will contain the total time in ns for each different ion on average
    diff_ions_number = [] # list that will contain the total number of different ions for each simulation on average

    for i in range(min, max): # loop through all the ions
        ion_df = pd.read_table(filename+str(i)+'-' +residue +'.dat', sep = r'\s+')# read in the distancedata for one ion, 
        # make it into a dataframe
        ion_df = ion_df.drop(['#Frame'], axis=1)
        files_list.append(ion_df.rename(columns={"Dis_00001": str(i)})) # add it to the list of dataframes

    ion_distances = pd.concat(files_list, axis=1) # make a new dataframe variable, 'ion_distances' that is the combined dataframe for all ions

    less_than_4 = ion_distances < 4 # make new data frame and tell whether ion_distances are less than 4 

    if statusfile == "none":
        contact_by_frame = less_than_4.sum(axis = 1) > 0 # take what 'less_than_4' spit out and sum through each row to see if there are any true values
        ions_by_frame = less_than_4.sum(axis = 1) # total number of ions in each frame
        
        #cut up columns into each independent simulation
        sim_len = contact_by_frame.count()/nsim # divide total number of frames by total number of simulations to get length of one simulation 
        for i in range(0, nsim): # loop through each independent simulation
            sim_contacts = contact_by_frame[int(sim_len*i):int(sim_len*(i+1))]
            avg_for_sim = np.mean(sim_contacts) # calculate fraction of simulation where there's an ion contact in each simulation
            sim_averages.append(avg_for_sim) # take simulation averages and combine into a list of fraction of simulation where there's an ion contact
            
            sim_ions = ions_by_frame[int(sim_len*i):int(sim_len*(i+1))]
            avg_ions_for_sim = np.mean(sim_ions) # calculate fraction of simulation where there's an ion contact in each simulation
            sim_ions_averages.append(avg_ions_for_sim) # take simulation averages and combine into a list of fraction of simulation where there's an ion contact
            
            # calculate the frac of time for each distinct ion
            sim_less_than_4 = less_than_4[int(sim_len*i):int(sim_len*(i+1))]
            sim_different_ions = sim_less_than_4.mean(axis = 0) # frac for each ion
            diff_ions_averages.append(sim_different_ions.mean()*(sim_len/100)) # mean time in ns for individual ion in this sim
            diff_ions_number.append((sim_different_ions > 0).sum()) # total number of ions in contact in this sim

    else:
        contact_by_frame = pd.DataFrame(columns=['contact', 'EC'])
        ions_by_frame = pd.DataFrame(columns=['ions', 'EC'])
        
        contact_by_frame['contact'] = less_than_4.sum(axis = 1) > 0 # take what 'less_than_4' spit out and sum through each row to see if there are any true values
        ions_by_frame['ions'] = less_than_4.sum(axis = 1) # total number of ions in each frame
        
        contact_by_frame['EC'] = (stat['status'] != 'unbound') & (stat['status'] != 'Bound') # tell whether frame is in encounter complex
        ions_by_frame['EC'] = (stat['status'] != 'unbound') & (stat['status'] != 'Bound') # tell whether frame is in encounter complex
        
        # add the EC column to the less_than_4 data dataframe to enable calculating frac for each distinct ion
        less_than_4['EC'] = contact_by_frame['EC']
        
        #cut up columns into each independent simulation
        sim_len = contact_by_frame['contact'].count()/nsim # divide total number of frames by total number of simulations to get length of one simulation 
    
        for i in range(0, nsim): # loop through each independent simulation
            sim_contacts = contact_by_frame.iloc[int(sim_len*i):int(sim_len*(i+1))]
            avg_for_sim = np.mean(sim_contacts[sim_contacts['EC']]['contact']) # calculate fraction of simulation where there's an ion contact in each simulation if in encounter complex
            sim_averages.append(avg_for_sim) # take simulation averages and combine into a list of fraction of simulation where there's an ion contact
    
            sim_ions = ions_by_frame.iloc[int(sim_len*i):int(sim_len*(i+1))]
            avg_ions_for_sim = np.mean(sim_ions[sim_ions['EC']]['ions']) # calculate fraction of simulation where there's an ion contact in each simulation
            sim_ions_averages.append(avg_ions_for_sim) # take simulation averages and combine into a list of fraction of simulation where there's an ion contact
            
            # calculate the frac of time for each distinct ion
            sim_less_than_4 = less_than_4.iloc[int(sim_len*i):int(sim_len*(i+1))]
            sim_different_ions = sim_less_than_4[sim_less_than_4['EC']].mean(axis = 0) # frac for each ion
            diff_ions_averages.append(sim_different_ions.mean()*(sim_len/100)) # mean time in ns for individual ion in this sim
            diff_ions_number.append((sim_different_ions > 0).sum()) # total number of ions in contact in this sim

    finalfile = [] # create new list 'finalfile'
    finalfile.append(np.mean(sim_averages)) # calculate the average fraction with an ion contact
    finalfile.append(stats.sem(sim_averages)) # take fraction of simulation where there's an ion contact data and calculate standard deviation
    
    finalfile.append(np.mean(sim_ions_averages)) # calculate the average fraction with an ion contact
    finalfile.append(stats.sem(sim_ions_averages)) # take fraction of simulation where there's an ion contact data and calculate standard deviation
    
    # save the avg time for each distinct ion
    finalfile.append(np.mean(diff_ions_averages)) # average time for an individual ion
    finalfile.append(stats.sem(diff_ions_averages)) # standard error
    # save the total number of distinct ions
    finalfile.append(np.mean(diff_ions_number)) # average number of individual ions
    finalfile.append(stats.sem(diff_ions_number)) # standard error
    
    return finalfile


# function to loop through all charged residues in AbpSH3 to find the fraction of the ensemble in contact with an ion
def all_ion_contacts(construct, path, protein_res_num, pos_ion_num, neg_ion_num, nsim, pos_ion, neg_ion, pos_first, res_names, statusfile):
    n_groups = len(res_names) # determine total number of residues
    NaCl_final_contact_fraction_df = pd.DataFrame(columns = ['res']) # create a dataframe for the 
                                                                    #final contact fractions for each residue
    NaCl_final_contact_fraction_df['res'] = res_names # put residue names into dataframe
    NaCl_final_contact_fraction_df = NaCl_final_contact_fraction_df.set_index('res') # convert residue names to the row indeces
    NaCl_avg_ions_df = NaCl_final_contact_fraction_df.copy() # create a dataframe for the average number of ions for each residue
    NaCl_avg_diff_ions_df = NaCl_final_contact_fraction_df.copy() # create a dataframe for the time for disinct ions for each residue
    NaCl_avg_num_diff_ions_df = NaCl_final_contact_fraction_df.copy() # create a dataframe for number of disinct ions for each residue
    
    for name in res_names: # loop through all the residue names
    
        # define parameters for our ion_contact function
        if name[0:3] == 'Lys' : # test if the residues is lysine
            filename = path + neg_ion + '_distances/lys_' + neg_ion + '_' # set the path and filename prefix for lysine
            residue = name[3:] # assign the residue number to the variable 'residue'
            if residue[0] == '(' : # test if there are parantheses in the residue number
                residue = residue[1:-1] # get rid of the parantheses if they are there
                if residue[0] == '-' : # test to see if there is a negative
                    residue = 'neg' + residue[1:] # change '-' to 'neg' in residue number
            if pos_first:
                min = protein_res_num + pos_ion_num + 1 # determine the first residue number for the negative ions
            else:
                min = protein_res_num + 1
            max = min + neg_ion_num # determine the last residue number for the negative ions
            
        if name[0:3] == 'Glu' : # test if the residue is glutamate
            filename = path + pos_ion + '_distances/glu_' + pos_ion + '_' # set path and prefix for glutamate
            residue = name[3:] # assign residue number
            if pos_first:
                min = protein_res_num + 1 # determine first positive ion res number
            else:
                min = protein_res_num + neg_ion_num + 1
            max = min + pos_ion_num # determine last positive ion res number
            
        if name[0:3] == 'Asp' : # test if the residue is aspartate
            filename = path + pos_ion + '_distances/asp_' + pos_ion + '_' # set path and prefix for aspartate
            residue = name[3:] # assign residue number
            if pos_first:
                min = protein_res_num + 1 # determine first positive ion res number
            else:
                min = protein_res_num + neg_ion_num + 1
            max = min + pos_ion_num # determine last positive ion res number
        
        # run the ion contact function for each residue
        NaCl_final_contact_fraction = ion_contact(filename, residue, min, max, nsim, statusfile) # run 'ion_contacts'
        NaCl_final_contact_fraction_df.loc[name, construct + ' avg'] = NaCl_final_contact_fraction[0] # store the fraction of ensemble calculated
        NaCl_final_contact_fraction_df.loc[name, construct + ' sem'] = NaCl_final_contact_fraction[1] # store the standard error calculated
        
        NaCl_avg_ions_df.loc[name, construct + ' avg'] = NaCl_final_contact_fraction[2] # store the fraction of ensemble calculated
        NaCl_avg_ions_df.loc[name, construct + ' sem'] = NaCl_final_contact_fraction[3] # store the standard error calculated

    
        NaCl_avg_diff_ions_df.loc[name, construct + ' avg'] = NaCl_final_contact_fraction[4] # store the avg time
        NaCl_avg_diff_ions_df.loc[name, construct + ' sem'] = NaCl_final_contact_fraction[5] # store the standard error
        
        NaCl_avg_num_diff_ions_df.loc[name, construct + ' avg'] = NaCl_final_contact_fraction[6] # store the number of ions
        NaCl_avg_num_diff_ions_df.loc[name, construct + ' sem'] = NaCl_final_contact_fraction[7] # store the standard error
    
    return NaCl_final_contact_fraction_df, NaCl_avg_ions_df, NaCl_avg_diff_ions_df, NaCl_avg_num_diff_ions_df


# function to find the average total number of cation and anion contacts as well as the total number of ions in contact with a charged residue
def total_ion_contacts(final_contact_fraction_df, avg_ions_df, res_names):
    
    # create a dataframe for the average number of cations and anions
    # to get the uncertainty, first we will need the square of it
    total_contacts = pd.DataFrame(0, index = pd.Series(["cation", "anion"]), 
                                  columns = ["avg", "semsq"], dtype=np.float64)

    for name in res_names: # loop through all the residue names
        
        if name[0:3] == 'Lys' : # test if the residues is lysine
            total_contacts.loc["anion","avg"] = total_contacts.loc["anion","avg"] + final_contact_fraction_df.loc[name].iloc[0]
            total_contacts.loc["anion","semsq"] = (total_contacts.loc["anion","semsq"] + 
                                                   np.square(final_contact_fraction_df.loc[name].iloc[1]))
            
        if (name[0:3] == 'Glu') or (name[0:3] == 'Asp') : # test if the residues is glu or asp
            total_contacts.loc["cation","avg"] = total_contacts.loc["cation","avg"] + final_contact_fraction_df.loc[name].iloc[0]
            total_contacts.loc["cation","semsq"] = (total_contacts.loc["cation","semsq"] + 
                                                    np.square(final_contact_fraction_df.loc[name].iloc[1]))
    
    total_contacts["sem"] = np.sqrt(total_contacts["semsq"])

    # find the average total number of cations and anions
    total_ions = pd.DataFrame(0, index = pd.Series(["cation", "anion"]), 
                                  columns = ["avg", "semsq"], dtype=np.float64)

    for name in res_names: # loop through all the residue names
        
        if name[0:3] == 'Lys' : # test if the residues is lysine
            total_ions.loc["anion","avg"] = total_ions.loc["anion","avg"] + avg_ions_df.loc[name].iloc[0]
            total_ions.loc["anion","semsq"] = (total_ions.loc["anion","semsq"] + 
                                                   np.square(avg_ions_df.loc[name].iloc[1]))
            
        if (name[0:3] == 'Glu') or (name[0:3] == 'Asp') : # test if the residues is glu or asp
            total_ions.loc["cation","avg"] = total_ions.loc["cation","avg"] + avg_ions_df.loc[name].iloc[0]
            total_ions.loc["cation","semsq"] = (total_ions.loc["cation","semsq"] + 
                                                    np.square(avg_ions_df.loc[name].iloc[1]))
    
    total_ions["sem"] = np.sqrt(total_ions["semsq"])

    
    return total_contacts, total_ions


# A function to run all ion contacts for multiple types of simulations
def compare_sims(comparison_sims, sim_colors, save_path, fig_path, label):

    full_ion_contact_data = pd.DataFrame()
    full_avg_ion_data = pd.DataFrame()
    full_avg_diff_ion_data = pd.DataFrame()
    full_avg_num_diff_ion_data = pd.DataFrame()
    
    full_total_contacts_data = pd.DataFrame()
    full_total_ions_data = pd.DataFrame()
    
    for construct in comparison_sims: # loop through the types of sims we want to compare
        # this part I should change to reading in from a file
        path = '/data/acarhart/salt_project/abp1_NaCl/analysis/'
        protein_res_num = 58 # number of residues in the protein
        pos_ion_num = 51 # number of positively charged ions in the simulation
        neg_ion_num = 39 # number of negatively charged ions in the simulation
        nsim = 10 # number of independent simulations
        pos_ion = 'Na'
        neg_ion = 'Cl'
        pos_first = True # postive ions are first in the topology file
        res_names = ('Glu7', 'Asp9', 'Asp11', 'Glu14', 'Asp15', 'Glu17', 'Glu22', 'Asp24', 'Lys25', 'Glu30', 'Asp33', 
                 'Asp34', 'Asp35', 'Glu40', 'Glu42', 'Lys43', 'Asp44', 'Lys47')
        
        statusfile = 'none'
        if construct == 'ArkA17_bound':
            path = '/data/kball/SH3_binding/salt_paper_analysis/AbpSH3_ArkA17_NaCl/'
            protein_res_num = 76 # number of residues in the protein
            pos_ion_num = 68 # number of positively charged ions in the simulation
            neg_ion_num = 61 # number of negatively charged ions in the simulation
            nsim = 10 # number of independent simulations
            pos_ion = 'Na'
            neg_ion = 'Cl'
            pos_first = True # negative ions are first in the topology file
            res_names = ('Glu7', 'Asp9', 'Asp11', 'Glu14', 'Asp15', 'Glu17', 'Glu22', 'Asp24', 'Lys25', 'Glu30', 'Asp33', 
                         'Asp34', 'Asp35', 'Glu40', 'Glu42', 'Lys43', 'Asp44', 'Lys47', 'Lys(6)', 'Lys(5)', 'Lys(3)', 
                         'Lys(-3)', 'Lys(-8)', 'Lys(-10)') # list all the charged residues
        elif construct == 'ArkA12_bound':
            path = '/data/kball/SH3_binding/salt_paper_analysis/AbpSH3_ArkA12_NaCl/'
            protein_res_num = 72 # number of residues in the protein
            pos_ion_num = 52 # number of positively charged ions in the simulation
            neg_ion_num = 43 # number of negatively charged ions in the simulation
            nsim = 10 # number of independent simulations
            pos_ion = 'Na'
            neg_ion = 'Cl'
            pos_first = False # negative ions are first in the topology file
            res_names = ('Glu7', 'Asp9', 'Asp11', 'Glu14', 'Asp15', 'Glu17', 'Glu22', 'Asp24', 'Lys25', 'Glu30', 'Asp33', 
                         'Asp34', 'Asp35', 'Glu40', 'Glu42', 'Lys43', 'Asp44', 'Lys47', 'Lys(3)', 
                         'Lys(-3)', 'Lys(-8)') # list all the charged residues
        elif construct == 'ArkA12_encounter':
            path = '/data/kball/SH3_binding/salt_paper_analysis/ArkA12_binding_NaCl/'
            protein_res_num = 72 # number of residues in the protein
            pos_ion_num = 193 # number of positively charged ions in the simulation
            neg_ion_num = 184 # number of negatively charged ions in the simulation
            nsim = 50 # number of independent simulations
            pos_ion = 'Na'
            neg_ion = 'Cl'
            pos_first = False # negative ions are first in the topology file
            res_names = ('Glu7', 'Asp9', 'Asp11', 'Glu14', 'Asp15', 'Glu17', 'Glu22', 'Asp24', 'Lys25', 'Glu30', 'Asp33', 
                         'Asp34', 'Asp35', 'Glu40', 'Glu42', 'Lys43', 'Asp44', 'Lys47', 'Lys(3)', 
                         'Lys(-3)', 'Lys(-8)') # list all the charged residues
            statusfile = '/data/kball/SH3_binding/salt_paper_analysis/ArkA12_binding_NaCl/ArkA12_status.csv'

        final_contact_fraction_df, avg_ions_df, avg_diff_ions_df, avg_num_diff_ions_df = all_ion_contacts(construct, path, protein_res_num, pos_ion_num, neg_ion_num, nsim, pos_ion, neg_ion, pos_first, res_names, statusfile)
        total_contacts, total_ions = total_ion_contacts(final_contact_fraction_df, avg_ions_df, res_names)
        
        # merge data frames
        if full_ion_contact_data.empty:
            full_ion_contact_data = final_contact_fraction_df
        else:
            full_ion_contact_data = full_ion_contact_data.join(final_contact_fraction_df)
            
        if full_avg_ion_data.empty:
            full_avg_ion_data = avg_ions_df
        else:
            full_avg_ion_data = full_avg_ion_data.join(avg_ions_df)
            
        if full_avg_diff_ion_data.empty:
            full_avg_diff_ion_data = avg_diff_ions_df
        else:
            full_avg_diff_ion_data = full_avg_diff_ion_data.join(avg_diff_ions_df)
            
        if full_avg_num_diff_ion_data.empty:
            full_avg_num_diff_ion_data = avg_num_diff_ions_df
        else:
            full_avg_num_diff_ion_data = full_avg_num_diff_ion_data.join(avg_num_diff_ions_df)
        
        if full_total_contacts_data.empty:
            full_total_contacts_data = total_contacts
        else:
            full_total_contacts_data = full_total_contacts_data.join(total_contacts, rsuffix=construct)    
            
        if full_total_ions_data.empty:
            full_total_ions_data = total_ions
        else:
            full_total_ions_data = full_total_ions_data.join(total_ions, rsuffix=construct)
    
    
    res_names = list(full_ion_contact_data.index.values)
    full_ion_contact_data = full_ion_contact_data.fillna(0)
    full_avg_ion_data = full_avg_ion_data.fillna(0)
    full_avg_diff_ion_data = full_avg_diff_ion_data.fillna(0)
    full_avg_num_diff_ion_data = full_avg_num_diff_ion_data.fillna(0)

    
    # write the data to files
    full_ion_contact_data.to_csv(save_path + 'full_ion_contact_data_' + label + '.csv')
    full_avg_ion_data.to_csv(save_path + 'full_avg_ion_data_' + label + '.csv')
    full_avg_diff_ion_data.to_csv(save_path + 'full_avg_time_diff_ion_data_' + label + '.csv')
    full_avg_num_diff_ion_data.to_csv(save_path + 'full_avg_num_diff_ion_data_' + label + '.csv')
    full_total_contacts_data.to_csv(save_path + 'full_total_contacts_data_' + label + '.csv')
    full_total_ions_data.to_csv(save_path + 'full_total_ions_data_' + label + '.csv')

    # Make the plot

    plt.style.use('default')
    plt.rcParams['figure.figsize'] = (1 , .5)
    plt.rc('axes', titlesize=7)     # fontsize of the axes title
    plt.rc('axes', labelsize=7)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=7)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=7)    # fontsize of the tick labels
    plt.rc('legend', fontsize=7)    # fontsize of the legend
    
    bar_width = 0.28
    opacity = 0.8
    plt.figure(figsize=(11,4)) # set size of figure, width and height
    
    x_values = np.arange(len(res_names))
    
    for idx, construct in enumerate(comparison_sims): # loop through the types of sims we want to compare
        rects = plt.bar(x_values+(idx*bar_width), full_ion_contact_data[construct + " avg"], bar_width, yerr = full_ion_contact_data[construct + " sem"], error_kw={'elinewidth': 1, "capsize": 1.5, "capthick": 0.5}, alpha=opacity, color=sim_colors[idx], label=construct.replace("_", " "))       
        
    plt.xlabel('Residue',fontsize=7)
    plt.ylabel('Fraction of Simulation Time', fontsize=7)
    #plt.title('Contact Between Charged Residues and Ions')
    plt.xticks(x_values + bar_width, res_names) # put the residue names on the x-axis for each bar
    
    plt.legend(fontsize=7,frameon=False,loc='upper right') # add a legend
    plt.xticks(rotation=30) # put the x-axis labels at an angle
    
    # Get the current figure
    fig = plt.gcf()
    
    # Get the current dimensions
    current_width, current_height = fig.get_size_inches()
    print(f"Current dimensions: {current_width} x {current_height} inches")
    
    # Define the new width
    new_width = 7
    
    # Calculate the new height to maintain the aspect ratio
    aspect_ratio = current_height / current_width
    new_height = new_width * aspect_ratio
    print(f"New dimensions: {new_width} x {new_height} inches")
    
    # Set the new figure size
    fig.set_size_inches(new_width, new_height)
    
    # save the figure to a file
    plt.savefig(fig_path, bbox_inches = 'tight', dpi = 1000)

    return full_ion_contact_data, full_avg_ion_data, full_total_contacts_data, full_total_ions_data


if __name__ =='__main__': 
    plt.switch_backend('agg')
    # colors: Apo = 'darkorange'
    # ArkA12_bound = '#332288'
    # ArkA17_bound = '#44AA99'
    # ArkA12_encounter = '#88CCEE'

    comparison_sims = ['Apo', 'ArkA12_encounter', 'ArkA12_bound']
    sim_colors = ['darkorange', '#88CCEE', '#332288']
    save_path = '/data/mcohen3/SH3/notebooks/'
    label = 'apo_12EC_12bound'
    fig_path = '/data/mcohen3/SH3/notebooks/800mM_NaCl_ion_contacts_' + label + '.png'

    compare_sims(comparison_sims, sim_colors, save_path, fig_path, label)

