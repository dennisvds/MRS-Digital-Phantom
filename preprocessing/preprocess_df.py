####################################################################################################
#                                       preprocess_df.py                                           #
####################################################################################################
#                                                                                                  #
# Authors: D.M.J. van de Sande (d.m.j.v.d.sande@tue.nl)                                            #
#                                                                                                  #
# Created: 26/01/2024                                                                              #
#                                                                                                  #
# Purpose: Script to generate and prepare a metabolite dataframe as input information for the      #
#          MRS/MRSI Digital Phantom. It take the MRS database from:                                #
#                                                                                                  #                       
#           Gudmundson AT, Koo A, Virovka A, et al.                                                #
#           Meta-analysis and open-source database for in vivo brain Magnetic Resonance            #
#           spectroscopy in health and disease.                                                    #
#           Analytical Biochemistry. 2023;676:115227. doi:10.1016/j.ab.2023.115227                 #          
#                                                                                                  #
#                                                                                                  #
####################################################################################################


####################################################################################################
#                                       Import packages                                            #
####################################################################################################
import pandas as pd
import numpy as np
import re

# Own
from utils.definitions import METABS, LABELS, METABS_TRANSLATION, CSF_DATA
from utils.supp_functions import fixed_effect, random_effect


def print_stats(df_old, df_new, step_name):
    """
    Prints statistics about the dataframes before and after a specific step.

    Parameters:
        df_old (pandas.DataFrame): The original dataframe.
        df_new (pandas.DataFrame): The modified dataframe.
        step_name (str): The name of the step being performed.

    Returns:
        None
    """
    print(f'''{step_name}... 
          Before: Unique articles - {len(df_old.Reference.unique())}, Entries - {len(df_old)} 
          After : Unique articles - {len(df_new.Reference.unique())}, Entries - {len(df_new)}
        ''')
    
def split_combined_metabs(df, combined_metab, metabs, ratios_WM, ratios_GM):
    """
    Splits a combined metabolite into individual metabolites based on white matter (WM) and gray matter (GM) ratios.

    Args:
        df (pandas.DataFrame): The input DataFrame containing metabolite data.
        combined_metab (str): The name of the combined metabolite to be split.
        metabs (list): A list of individual metabolite names to be created.
        ratios_WM (list): A list of ratios for splitting the combined metabolite in white matter.
        ratios_GM (list): A list of ratios for splitting the combined metabolite in gray matter.

    Returns:
        pandas.DataFrame: The updated DataFrame with the combined metabolite split into individual metabolites.

    """
    print(f'Splitting {combined_metab} into {metabs}...')

    # select all entries with the combined metabolite
    df_combined = df[df['Metabolite'] == combined_metab]

    # split in white and gray matter
    df_combined_WM = df_combined[df_combined['Tissue'] == 'WM'].copy()
    df_combined_GM = df_combined[df_combined['Tissue'] == 'GM'].copy()

    result_dfs = []

    for metab, ratio_WM, ratio_GM in zip(metabs, ratios_WM, ratios_GM):
        # split concentrations according to ratios for WM
        df_combined_WM.loc[:, ['IU_u', 'mM_u', 'IU_std', 'mM_std']] *= ratio_WM
        df_metab_WM = df_combined_WM.copy()
        df_metab_WM['Metabolite'] = metab

        # split concentrations according to ratios for GM
        df_combined_GM.loc[:, ['IU_u', 'mM_u', 'IU_std', 'mM_std']] *= ratio_GM
        df_metab_GM = df_combined_GM.copy()
        df_metab_GM['Metabolite'] = metab

        result_dfs.extend([df_metab_WM, df_metab_GM])

    # replace combined metabolite entries with resulting metabolites
    df = df[df['Metabolite'] != combined_metab]
    df = pd.concat([df] + result_dfs, ignore_index=True)

    return df

def change_metab_names(df):
    """
    Change the names of metabolites in the given DataFrame based on a predefined translation dictionary.

    Args:
        df (pandas.DataFrame): The DataFrame containing the metabolite data.

    Returns:
        pandas.DataFrame: The DataFrame with updated metabolite names.

    """
    for metab in df.Metabolite.unique():
        if metab not in METABS:
            print(f"Metabolite '{metab}' is not in METABS...")
            try:
                print(f"Change metabolite name '{metab}' to '{METABS_TRANSLATION[metab]}'...")
                df.loc[df.Metabolite == metab, 'Metabolite'] = METABS_TRANSLATION[metab]
            except KeyError:
                print(f"Metabolite '{metab}' is not in METABS_TRANSLATION...")
    return df.reset_index(drop=True)
    
def process_age(age):
    """
    Process the age value and return a float.

    Parameters:
    age (str or numeric): The age value to be processed.

    Returns:
    float: The processed age value.

    Raises:
    ValueError: If the age value cannot be converted to a float.

    """
    if isinstance(age, str):
        # If the age is a string, try to extract a range
        match = re.match(r'(\d+)-(\d+)', age)

        if match:
            # If it's a range, calculate the average
            lower_bound, upper_bound = map(int, match.groups())
            return (lower_bound + upper_bound) / 2
        else:
            # If it's not a range, try converting to float
            try:
                return float(age)
            except ValueError:
                # If conversion fails, raise an error
                raise ValueError(f"Unable to convert '{age}' to a float.")
    else:
        # If it's already a number, return the number
        return float(age)

def process_age_column(df):
    """
    Process the 'Age' column of a DataFrame by applying the process_age function to each value.

    Args:
        df (pandas.DataFrame): The DataFrame containing the 'Age' column.

    Returns:
        pandas.DataFrame: The DataFrame with the 'Age' column processed.

    """
    # Apply the process_age function to the 'Age' column
    df['Age'] = df['Age'].apply(process_age)
    return df

def preprocess_df(df, group, fraction_boundary=0.6, age_range=[18, 60]):
    """
    Preprocesses the given DataFrame by performing several data transformations and filtering operations.

    Args:
        df (pandas.DataFrame): The input DataFrame to be preprocessed.
        group (str or list): The group(s) to filter on.
        fraction_boundary (float, optional): The boundary value used to assign tissue labels based on GM/WM fractions. Defaults to 0.6.
        age_range (list, optional): The age range to select on. Defaults to [18, 60].

    Returns:
        pandas.DataFrame: The preprocessed DataFrame.

    """
    
    # Filter on a specific group
    print_stats(df, df[df.Group.isin(group)], f'Filtering on "{group}" group')
    df = df[df.Group.isin(group)]

    # Drop entries where no information about GM/WM is present
    columns_to_check = ['Tissue', 'GM', 'GM_Std', 'WM', 'WM_Std']
    print_stats(df, df.dropna(subset=columns_to_check, how='all'), 'Dropping entries where no information about GM/WM is present')
    df = df.dropna(subset=columns_to_check, how='all')

    # Change percentages to fractions
    print('Change percentages to fractions...')
    tissue_columns = ['GM', 'WM']
    std_columns = ['GM_Std', 'WM_Std']

    for idx, column in enumerate(tissue_columns):
        value_column = df[column]
        std_column = std_columns[idx]

        # Apply the transformation to the value column
        df[column] = value_column.apply(lambda x: x / 100 if x > 1 else x)

        # Apply the same transformation to the corresponding standard deviation column only if the value column is transformed
        if (value_column > 1).any():
            df[std_column] = df[std_column].apply(lambda x: x / 100 if x > 1 else x)

    # Assign tissue value based on GM/WM fractions (leaving already assigned tissue labels untouched)
    print(f'Assign tissue value based on GM/WM fractions (boundary: {fraction_boundary})...')
    df['Tissue'] = df.apply(lambda x: 'GM' if x.GM >= fraction_boundary else ('WM' if x.WM >= fraction_boundary else x.Tissue), axis=1)

    # Drop all entries that are not WM or GM
    print_stats(df, df[(df.Tissue == 'WM') | (df.Tissue == 'GM')], 'Drop all entries that are not WM or GM')
    df = df[(df.Tissue == 'WM') | (df.Tissue == 'GM')]

    # Change NaN values to 0
    print('Change NaN values to 0...')
    df = df.fillna(0)

    # Select on age range
    df = process_age_column(df)
    print_stats(df, df[(df.Age >= age_range[0]) & (df.Age <= age_range[1])], f'Select age range of {age_range[0]}-{age_range[1]} years')
    df = df[(df.Age >= age_range[0]) & (df.Age <= age_range[1])]

    # Remove all MM entries
    print_stats(df, df[~df.Name.str.contains('MM')], 'Remove all MM entries')
    df = df[~df.Name.str.contains('MM')]

    # Change 'Name' to 'metabolite'
    print('Change "Name" to "Metabolite"...')
    df = df.rename(columns={'Name': 'Metabolite'})

    return df.reset_index(drop=True)



def load_mrs_database(groups=['Healthy', 'Control'], fraction_boundary=0.6, age_range=[18, 60], tesla=3.0):
    """
    Load and preprocess the MRS database.

    Parameters:
    - groups (list): A list of group names to filter the database by. Default is ['Healthy', 'Control'].

    Returns:
    - dfConcs (DataFrame): Preprocessed concentration dataframe.
    - dfT2 (DataFrame): Preprocessed T2 dataframe.
    """

    print('LOADING MRS DATABASE...')
    # Load MRS Database
    dfCRef = pd.read_excel('data/metabolites/MRS_Database.xlsx', 'Refs_Concs')
    dfCVal = pd.read_excel('data/metabolites/MRS_Database.xlsx', 'Values_Concs')

    dfT2Ref = pd.read_excel('data/metabolites/MRS_Database.xlsx', 'Refs_T2')
    dfT2Val = pd.read_excel('data/metabolites/MRS_Database.xlsx', 'Values_T2')

    # Merge information from both sheets into one dataframe
    dfConcs = pd.merge(dfCVal[['Name', 'Reference', 'ID', 'DOI', 'tCr_u', 'tCr_std', 'IU_u', 'IU_std', 'mM_u', 'mM_std']], 
                        dfCRef[['Reference', 'ID', 'Group', 'Tissue', 'Title', 'GM', 'GM_Std', 'WM', 'WM_Std', 'N_Total', 'Female', 'Male', 'Age', 'Age_Std']], 
                        on=['Reference', 'ID'], how='inner').reset_index(drop=True)
    
    dfT2 = pd.merge(dfT2Val[['Name', 'Reference', 'ID', 'DOI', 'Tesla', 'T2', 'StdDev']], 
                    dfT2Ref[['Reference', 'ID', 'Group', 'Tissue', 'Title', 'GM', 'GM_Std', 'WM', 'WM_Std', 'N_Total', 'Female', 'Male', 'Age', 'Age_Std', 'Subject']],
                    on=['Reference', 'ID'], how='inner').reset_index(drop=True)
    
    # Preprocess concentration dataframe
    print('Preprocess concentration dataframe...')
    dfConcs = preprocess_df(dfConcs, groups, fraction_boundary, age_range)

    # Preprocess T2 dataframe
    print('Preprocess T2 dataframe...')
    dfT2 = dfT2[dfT2.Subject == 'Human'].drop(columns=['Subject']).reset_index(drop=True)
    dfT2 = preprocess_df(dfT2, groups, fraction_boundary, age_range)

    # Filter on magnetic field strength
    print('Filter T2 data on magnetic field strength...')
    print_stats(dfT2, dfT2[dfT2.Tesla == tesla].reset_index(drop=True), f'Filtering on {tesla} Tesla')
    dfT2 = dfT2[dfT2.Tesla == tesla].reset_index(drop=True)

    # Split combined metabolites
    # Based on:
    # [tNAA]               --> Pouwels PJW, Frahm J. Differential distribution of NAA and NAAG in human brain as determined by quantitative localized proton MRS. NMR in Biomedicine. 
    #                          1997;10(2):73-78. doi:10.1002/(SICI)1099-1492(199704)10:2<73::AID-NBM448>3.0.CO;2-4
    # [Glx], [tCho], [tCr] --> Govindaraju V, Young K, Maudsley AA. Proton NMR chemical shifts and coupling constants for brain metabolites. NMR in Biomedicine. 
    #                          2000;13(3):129-153. doi:10.1002/1099-1492(200005)13:3<129::AID-NBM619>3.0.CO;2-V

    dfConcs = split_combined_metabs(dfConcs, combined_metab='tCr', metabs=['Cr', 'PCr'], ratios_WM=[0.5, 0.5], ratios_GM=[0.5, 0.5])
    dfConcs = split_combined_metabs(dfConcs, combined_metab='tNAA', metabs=['NAA', 'NAAG'], ratios_WM=[0.75, 0.25], ratios_GM=[0.9, 0.1])
    dfConcs = split_combined_metabs(dfConcs, combined_metab='Glx', metabs=['Glu', 'Gln'], ratios_WM=[0.66, 0.33], ratios_GM=[0.66, 0.33])
    dfConcs = split_combined_metabs(dfConcs, combined_metab='tCho', metabs=['PCh', 'GPC'], ratios_WM=[0.375, 0.625], ratios_GM=[0.375, 0.625])

    # Change metabolite names to definition names
    print('Change metabolite names in concentrations dataframe to names in definitions.py...')
    dfConcs = change_metab_names(dfConcs)
    print('Change metabolite names in T2 dataframe to names in definitions.py...')
    dfT2 = change_metab_names(dfT2)

    # Collapse mM and IU concentrations
    print('Collapse mM and IU concentrations...')
    dfConcs.IU_u  = dfConcs.IU_u.values + dfConcs.mM_u.values
    dfConcs.IU_std  = dfConcs.IU_std.values + dfConcs.mM_std.values
    
    print('DONE LOADING MRS DATABASE!')
    print(f'''Properties of dfConcs:
        Shape                       : {dfConcs.shape}
        Number of unique references : {len(dfConcs.Reference.unique())}
        Number of unique metabolites: {len(dfConcs.Metabolite.unique())}
        Metabolites with WM info    : {sorted(dfConcs[dfConcs.Tissue == 'WM'].Metabolite.unique())}
        Metabolites with GM info    : {sorted(dfConcs[dfConcs.Tissue == 'GM'].Metabolite.unique())}
    ''')

    print(f'''Properties of dfT2:
          Shape                       : {dfT2.shape}
          Number of unique references : {len(dfT2.Reference.unique())}
          Number of unique metabolites: {len(dfT2.Metabolite.unique())}
          Metabolites with WM info    : {sorted(dfT2[dfT2.Tissue == 'WM'].Metabolite.unique())}
          Metabolites with GM info    : {sorted(dfT2[dfT2.Tissue == 'GM'].Metabolite.unique())}
                                         ''')

    return dfConcs, dfT2

def calculate_metab_conc(df, metab, tissue):
    """
    Calculate the mean concentration and standard deviation of a metabolite in a specific tissue.

    Parameters:
    - df (pandas.DataFrame): The input dataframe containing metabolite data.
    - metab (str): The name of the metabolite.
    - tissue (str): The name of the tissue.

    Returns:
    - conc_mean (float): The mean concentration of the metabolite in the tissue.
    - conc_std (float): The standard deviation of the metabolite concentration in the tissue.
    """

    # Calculate mean concentration
    df_metab = df[df.Metabolite == metab].reset_index(drop=True)
    df_metab = df_metab[df_metab.Tissue == tissue].reset_index(drop=True)

    # Remove entries where no IU_u and IU_std are 0
    df_metab = df_metab[(df_metab.IU_u != 0) & (df_metab.IU_std != 0)].reset_index(drop=True)

    # Calculate weight column
    df_metab['Weight'] = 1 / df_metab.IU_std**2

    # Calculate random and fixed effects
    # Calculate mean concentration
    if len(df_metab) == 0:
        conc_mean, conc_std = np.nan, np.nan
    elif len(df_metab) == 1:
        conc_mean = df_metab['IU_u'].values[0]
        conc_std = df_metab['IU_std'].values[0]
    else:
        random_effects = random_effect(df_metab, Tcol='IU_u', Wcol='Weight')
        conc_mean = random_effects['M']
        conc_std = random_effects['SE']

    return conc_mean, conc_std

def calculate_metab_relaxation(df, metab, tissue, tesla=3.0):
    """
    Calculate the mean relaxation time for a specific metabolite and tissue at a given tesla value.

    Parameters:
    - df (pandas.DataFrame): The input DataFrame containing relaxation time data.
    - metab (str): The name of the metabolite.
    - tissue (str): The name of the tissue.
    - tesla (float, optional): The tesla value. Default is 3.0.

    Returns:
    - t2_mean (float): The mean relaxation time for the specified metabolite, tissue, and tesla value.
    """

    # Calculate mean relaxation time
    df_metab = df[df.Tesla == tesla].reset_index(drop=True)
    df_metab = df_metab[df_metab.Metabolite == metab].reset_index(drop=True)
    df_metab = df_metab[df_metab.Tissue == tissue].reset_index(drop=True)

    # Calculate weight column
    df_metab['Weight'] = 1 / df_metab.StdDev**2

    # Calculate random and fixed effects
    if len(df_metab) == 0:
        t2_mean, t2_std = np.nan, np.nan
    elif len(df_metab) == 1:
        t2_mean = df_metab['T2'].values[0]
        t2_std = df_metab['StdDev'].values[0]
    else:
        random_effects = random_effect(df_metab, Tcol='T2', Wcol='Weight')
        t2_mean = random_effects['M']
        t2_std = random_effects['SE']
        
    return t2_mean


# Create a list of dictionaries with random data for each metabolite, label, and label name
def create_metab_df(labels=['Background', 'WM', 'GM', 'CSF'], groups=['Healthy', 'Control'], 
                    fraction_boundary=0.6, age_range=[18, 60], tesla=3.0, save=False):
    """
    Create a metabolite dataframe based on concentration and relaxation data.

    Args:
        labels (list, optional): List of tissue labels. Defaults to ['Background', 'WM', 'GM', 'CSF'].
        save (bool, optional): Flag to save the metabolite dataframe. Defaults to False.

    Returns:
        tuple: A tuple containing the metabolite dataframe, concentration dataframe, and relaxation dataframe.
    """

    concs_df, t2_df = load_mrs_database(groups, fraction_boundary, age_range, tesla)
    # Unique metabolite names in both dataframes
    unique_metabs = np.concatenate([concs_df.Metabolite.unique(), t2_df.Metabolite.unique()])
    unique_metabs = np.unique(unique_metabs)
    # Remove MMs from unique_metabs
    unique_metabs = [metab for metab in unique_metabs if 'MM' not in metab]

    data = []
    for metab in unique_metabs:
        for label_id, label in enumerate(labels):
            # Background label
            if label_id == 0:
                conc_mean, conc_std, t2, t1 = 0, 0, 0, 0

            # GM & WM labels
            if label_id==1 or label_id==2:
                conc_mean, conc_std = calculate_metab_conc(concs_df, metab, label)
                t2 = calculate_metab_relaxation(t2_df, metab, label)
                t1 = np.nan  #TODO: add T1 information in the database

            # CSF label
            if label_id==3:
                if metab in CSF_DATA.keys():
                    conc_mean, conc_std, t1, t2 = CSF_DATA[metab]
                else:
                    conc_mean, conc_std = np.nan, np.nan
                    t2 = np.nan
                    t1 = np.nan

            entry = {
                'Metabolite': metab,
                'Label': label_id,
                'Tissue': label,
                'Conc_mean': np.round(conc_mean, 2),
                'Conc_std': np.round(conc_std, 2),
                'T1': np.round(t1, 2),
                'T2': np.round(t2, 2),
                }
                
            data.append(entry)

    # Create a pandas DataFrame from the list of dictionaries
    metab_df = pd.DataFrame(data)
    # Sort by metabolite and label
    metab_df = metab_df.sort_values(by=['Metabolite', 'Label']).reset_index(drop=True)
    metab_df = metab_df.sort_values(by='Metabolite', key=lambda x: x.str.lower()).reset_index(drop=True)
    
    if save:
        # Save metabolite df
        metab_df.to_csv('data/metabolites/metab_df.csv', index=False)
        concs_df.to_csv('data/metabolites/concs_df.csv', index=False)
        t2_df.to_csv('data/metabolites/t2_df.csv', index=False)
        print('Saved metabolite dataframe!')

    return metab_df, concs_df, t2_df


if __name__ == '__main__':
    # Calculate metabolite dataframe with mean concentrations
    metab_df, concs_df, t2_df = create_metab_df(labels=['Background', 'WM', 'GM', 'CSF'], groups=['Healthy', 'Control'], 
                                                fraction_boundary=0.6, age_range=[18, 60], save=False)
