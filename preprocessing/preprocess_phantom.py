####################################################################################################
#                                    preprocess_phantom.py                                         #
####################################################################################################
#                                                                                                  #
# Authors: D.M.J. van de Sande (d.m.j.v.d.sande@tue.nl)                                            #
#                                                                                                  #
# Created: 26/01/2024                                                                              #
#                                                                                                  #
# Purpose: Script to generate and preprocess the metabolite dataframe for the Digital Phantom.     #
#                                                                                                  #
####################################################################################################


####################################################################################################
#                                       Import packages                                            #
####################################################################################################
from tqdm import tqdm
import numpy as np
import nibabel as nib
import json
import os


def create_phantom_data(subject, metab_df, path2phantom, save=False):
    """
    Creates phantom data based on the given subject and metabolite dataframe.
    Args:
        subject (str): The subject for which the phantom data is created.
        metab_df (pandas.DataFrame): The metabolite dataframe containing the data.
        path2phantom (str): The path to save the phantom data.
        save (bool, optional): Whether to save the phantom data. Defaults to False.
    Returns:
        tuple: A tuple containing the phantom array and the metabolite name to ID mapping.
    """

    print('Creating phantom data...')
    phantom_array, metabolite_name_to_id = create_phantom_array(subject, metab_df)

    if save:
        print('Saving phantom data...')
        save_phantom_nii_files(path2phantom, phantom_array, metabolite_name_to_id, subject)
        print(f'Phantom data saved successfully to "{path2phantom}"')
    
    return phantom_array, metabolite_name_to_id


def create_phantom_array(subject, metab_df):
    """
    Create a phantom array based on the subject and metabolite dataframe.
    Args:
        subject (dict): A dictionary containing the subject data.
        metab_df (pandas.DataFrame): A dataframe containing metabolite information.
    Returns:
        tuple: A tuple containing the phantom array and the metabolite name to ID mapping.
    """ 
    # Extract label data from the subject dictionary
    label_data = subject['labels'].data[0]

    # Initialize the phantom dictionary and a mapping for metabolite names to numeric IDs
    phantom_dict = {}
    metabolite_name_to_id = {}

    print('Creating phantom dictionary...')
    for index, row in tqdm(metab_df.iterrows(), total=len(metab_df)):
        # Extract relevant information from each row of the dataframe
        label = row['Label']
        metabolite_name = row['Metabolite']
        concentration_mean = row['Conc_mean']
        concentration_std = row['Conc_std']
        t1_value = row['T1']
        t2_value = row['T2']

        # Find voxel indices where the label matches the current label in the row
        voxel_indices = np.where(label_data == label)

        # Assign a unique numeric ID to the metabolite if it hasn't been assigned already
        if metabolite_name not in metabolite_name_to_id:
            metabolite_name_to_id[metabolite_name] = len(metabolite_name_to_id)

        # Iterate over all matching voxel indices
        for x, y, z in zip(voxel_indices[0], voxel_indices[1], voxel_indices[2]):
            # Initialize the list for this voxel in the phantom dictionary if not already done
            if (x, y, z) not in phantom_dict:
                phantom_dict[(x, y, z)] = []

            # Append the metabolite information to the list for this voxel
            metabolite_info = [metabolite_name_to_id[metabolite_name], concentration_mean, concentration_std, t1_value, t2_value]
            phantom_dict[(x, y, z)].append(metabolite_info)
    

    # Determine the dimensions of the 3D label data
    x_size, y_size, z_size = label_data.shape
    # Determine the maximum number of metabolites per voxel
    max_metabolites_per_voxel = max(len(metabolites) for metabolites in phantom_dict.values())
    # Find the length of the metabolite information array (currently 5: [metabolite_id, conc_mean, conc_std, T1, T2])
    metabolite_info_size = len(metabolite_info)

    # Initialize the 5D phantom array (x, y, z, max_metabolites_per_voxel, metabolite_info_size)
    phantom_array = np.zeros((x_size, y_size, z_size, max_metabolites_per_voxel, metabolite_info_size), dtype=float)

    # Populate the phantom array with metabolite information
    print('Populating phantom array...')
    for (x, y, z), metabolite_infos in tqdm(phantom_dict.items()):
        for i, metabolite_info in enumerate(metabolite_infos):
            phantom_array[x, y, z, i, :] = metabolite_info
    
    # Return the final 5D phantom array and the metabolite name to ID mapping
    return phantom_array, metabolite_name_to_id


def save_phantom_nii_files(path2phantom, phantom_array, metabolite_name_to_id, subject):
    """
    Saves phantom and subject data as NIfTI files, including metadata.

    Parameters:
    - path2phantom (str): The directory where the NIfTI files will be saved.
    - phantom_array (numpy.ndarray): A 5D array representing the phantom data.
    - metabolite_name_to_id (dict): A dictionary mapping metabolite names to their numeric identifiers.
    - subject (dict): A dictionary containing subject data with 'labels', 'T1', 'T2', etc.
    
    Returns:
    - None (only saves the NIfTI files to the user-defined path).
    """
    
    # Initialize lists to store NIfTI file objects and their corresponding filenames
    nii_filenames = []
    nii_images = []

    # Extract the affine transformation matrix from the subject's label data
    affine_matrix = subject['labels'].affine

    # Prepare metadata to include in the NIfTI file header
    metadata = {
        'dim_info': {
            'DIM_1': 'x',
            'DIM_2': 'y',
            'DIM_3': 'z',
            'DIM_4': 'metabolite_id',
            'DIM_5': 'metabolite_info'
        },
        'metab_dim_info': {
            0: 'metabolite_id',
            1: 'conc_mean',
            2: 'conc_std',
            3: 't1',
            4: 't2'
        },
        'metab_mapping': metabolite_name_to_id
    }

    # Convert the metadata dictionary to a JSON string format
    metadata_json = json.dumps(metadata)

    # Create a NIfTI image for the phantom array with the associated affine matrix
    phantom_nii_image = nib.Nifti1Image(phantom_array, affine_matrix)
    # Append the metadata JSON as a NIfTI header extension
    phantom_nii_image.header.extensions.append(nib.nifti1.Nifti1Extension('comment', metadata_json.encode()))
    # Store the NIfTI image and its corresponding filename
    nii_images.append(phantom_nii_image)
    nii_filenames.append('phantom_metab.nii.gz')

    # Loop through each entry in the subject dictionary to create NIfTI files
    for key, value in subject.items():
        # Generate a filename for each NIfTI file
        filename = f'phantom_{key}.nii.gz'
        # Create a NIfTI image for the subject data
        nii_image = nib.Nifti1Image(value.data[0].numpy(), value.affine)
        # Store the NIfTI image and its corresponding filename
        nii_images.append(nii_image)
        nii_filenames.append(filename)

    # Save all the generated NIfTI files to the specified directory
    print('Saving NIfTI files...')
    for nii_image, filename in tqdm(zip(nii_images, nii_filenames), total=len(nii_filenames)):
        file_path = os.path.join(path2phantom, filename)
        nib.save(nii_image, file_path)
