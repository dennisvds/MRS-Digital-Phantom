####################################################################################################
#                                       load_skeleton.py                                           #
####################################################################################################
#                                                                                                  #
# Authors: D.M.J. van de Sande (d.m.j.v.d.sande@tue.nl)                                            #  
#                                                                                                  #                                                                                               
# Created: 30/06/22                                                                                #
#                                                                                                  #
# Purpose: Load and process files for the anatomy part of the Digital MRS Phantom.                 #
#           Currently implemented 'BigBrainMR' and 'MRiLab' loading functions.                     #
#                                                                                                  #
####################################################################################################

#*************#
#   imports   #
#*************#
import os
import torch
import numpy as np
import torchio as tio

def load_skeleton_as_subject(path2skeleton, skeleton_name, resolution=None):

    # Load skeleton data based on the requested skeleton
    if skeleton_name == 'BigBrainMR':
        subject = load_bigbrain_files(path2skeleton, resolution)
    elif skeleton_name == 'MRiLab':
        subject = load_mrilab_files(path2skeleton, resolution)
    else:
        raise NotImplementedError('The requested skeleton is not implemented yet.')

    return subject

def load_bigbrain_files(path2skeleton, resolution=None):
        
    subject = tio.Subject(
            labels=tio.LabelMap(os.path.join(path2skeleton, 'BigBrainMR_Labels_400um.nii.gz')),
            T1=tio.ScalarImage(os.path.join(path2skeleton, 'BigBrainMR_T1map_400um.nii.gz')),
            R2s=tio.ScalarImage(os.path.join(path2skeleton, 'BigBrainMR_T2smodel_R2smap_400um.nii.gz')),
            T1wImage=tio.ScalarImage(os.path.join(path2skeleton, 'BigBrainMR_T1weighted_400um.nii.gz')),
        )
    
    # Resample to isotropic resolution if requested
    if resolution is not None:
        subject = tio.Resample(resolution)(subject)

    # Preprocess the labels: 
    # only use background (0), white matter (1), gray matter (2), and CSF (3) set rest to background
    labels = subject['labels'].data[0]
    labels[labels > 3] = 0
    labels = torch.Tensor(labels).unsqueeze(0)                          # Convert back to tensor

    # Preprocess T1 map
    t1_data = subject['T1'].data[0].numpy().astype(float)               # Convert to numpy array
    t1_data = t1_data * 1e-3                                            # Convert to seconds for BigBrainMR
    t1_data = scale_T1_map(t1_data, subject['labels'].data[0].numpy())  # Scale T1 map to 3T using Rooney et al. (2007)
    t1_data = torch.Tensor(t1_data).unsqueeze(0)                        # Convert back to tensor      
    
    # Preprocess T2s map
    r2s_data = subject['R2s'].data[0].numpy().astype(float)             # Convert to numpy array
    r2s_data_copy = r2s_data.copy()                                     # Create a copy of the original array
    r2s_data_copy= r2s_data_copy / 100                                  # Convert to 1/seconds
    t2s_data = np.zeros_like(r2s_data_copy)                             # Create an array for T2s
    non_zero_indices = np.where(r2s_data_copy != 0)                     # Find non-zero indices
    t2s_data[non_zero_indices] = 1 / r2s_data_copy[non_zero_indices]    # Calculate T2s
    t2s_data = scale_T2_map(t2s_data, tesla=3.0)                        # Scale T2s map to 3T
    t2s_data = torch.Tensor(t2s_data).unsqueeze(0)                      # Convert back to tensor

    # Create a new subject with the preprocessed data
    subject_new = tio.Subject(
        labels=tio.LabelMap(tensor=labels, affine=subject['labels'].affine),
        T1=tio.ScalarImage(tensor=t1_data, affine=subject['T1'].affine),
        T2s=tio.ScalarImage(tensor=t2s_data, affine=subject['R2s'].affine),
        T1wImage=subject['T1wImage']
    ) 

    return subject_new

def load_mrilab_files(path2skeleton, resolution=None):
    subject = tio.Subject(
            labels=tio.LabelMap(os.path.join(path2skeleton, 'labels.nii.gz')),
            T1=tio.ScalarImage(os.path.join(path2skeleton, 'T1_map.nii.gz')),
            T2s=tio.ScalarImage(os.path.join(path2skeleton, 'T2s_map.nii.gz')),
            T2=tio.ScalarImage(os.path.join(path2skeleton, 'T2_map.nii.gz')),
            PD=tio.ScalarImage(os.path.join(path2skeleton, 'Rho_map.nii.gz')),
        )
    
    # Resample to isotropic resolution if requested
    if resolution is not None:
        subject = tio.Resample(resolution)(subject)
    
    # Preprocess the labels: 
    # Using background (0), white matter (1), gray matter (2), CSF (3), fat (4), and skull (5) set rest to background
    labels = subject['labels'].data[0]
    old_labels = labels.clone()
    labels[old_labels == 1] = 5 
    labels[old_labels == 2] = 1
    labels[old_labels == 3] = 2
    labels[old_labels == 4] = 3
    labels[old_labels == 5] = 4
    labels[labels > 5] = 0
    subject['labels'].data[0] = labels

    return subject

def scale_T1_map(t1_map, labels, tesla=3.0):
    # Scaling factors from Rooney et al. (2007), doi:10.1002/mrm.21122
    # Scaling 7T maps to 'tesla' T maps (default: 3T)

    # Scaling of WM
    t1_wm = t1_map[labels == 1]
    t1_wm_new = t1_wm * (tesla**0.382 / 7.0**0.382)

    # Scaling of GM
    t1_gm = t1_map[labels == 2]
    t1_gm_new = t1_gm * (tesla**0.376 / 7.0**0.376)

    # Scaling of CSF 
    t1_csf = t1_map[labels == 3]
    t1_csf_new = t1_csf * 1.0

    t1_new = t1_map.copy()
    t1_new[labels == 1] = t1_wm_new
    t1_new[labels == 2] = t1_gm_new
    t1_new[labels == 3] = t1_csf_new

    return t1_new

def scale_T2_map(t2_map, tesla=3.0):
    # Scaling is done by assuming a linear relationship between 1/T2 and field strength
    # Scaling 7T maps to 'tesla' T maps (default: 3T)

    t2_new = t2_map * (7.0 / tesla)
    
    return t2_new
