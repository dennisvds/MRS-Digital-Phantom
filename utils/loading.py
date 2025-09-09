import os
import numpy as np
import torch
import torchio as tio

from utils.auxillary import scale_T1_map, scale_T2_map


def load_skeleton_as_subject(path2skeleton, skeleton_name):

    # Load skeleton data based on the requested skeleton
    if skeleton_name == 'BigBrainMR':
        subject, label_names = load_bigbrain_files(path2skeleton)
    elif skeleton_name == 'MRiLab':
        subject, label_names = load_mrilab_files(path2skeleton)
    elif skeleton_name == 'Custom':
        subject, label_names = load_custom_files(path2skeleton)
    else:
        raise NotImplementedError('The requested skeleton is not implemented yet.')

    return subject, label_names

def load_bigbrain_files(path2skeleton):
        
    subject = tio.Subject(
            labels=tio.LabelMap(os.path.join(path2skeleton, 'BigBrainMR_Labels_400um.nii.gz')),
            T1=tio.ScalarImage(os.path.join(path2skeleton, 'BigBrainMR_T1map_400um.nii.gz')),
            R2s=tio.ScalarImage(os.path.join(path2skeleton, 'BigBrainMR_T2smodel_R2smap_400um.nii.gz')),
            T1wImage=tio.ScalarImage(os.path.join(path2skeleton, 'BigBrainMR_T1weighted_400um.nii.gz')),
        )
    
    # Create whole brain mask
    labels = subject['labels'].data[0]
    mask = (labels > 0).float()
    mask = mask.unsqueeze(0)  # Add a channel dimension

    # Original 20-label map from subject
    labels = subject['labels'].data[0]  # shape: (H, W, D)
    # Convert to numpy array for easier manipulation
    labels_np = labels.cpu().numpy()
    # Create a new array for the simplified labels
    simplified_labels = np.zeros_like(labels_np)
    # Assign new simplified labels
    # White Matter (1): original labels 1 and 4
    simplified_labels[np.isin(labels_np, [1, 4])] = 1
    # Gray Matter (2): original labels 2, 5–20
    simplified_labels[np.isin(labels_np, list(range(2, 21)))] = 2
    # CSF (3): original label 3
    simplified_labels[labels_np == 3] = 3
    # Convert back to tensor
    labels = torch.tensor(simplified_labels).unsqueeze(0).float()  # shape: (1, H, W, D)
    label_names = ['Background', 'WM', 'GM', 'CSF']

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
        T1wImage=subject['T1wImage'],
        brain_mask=tio.LabelMap(tensor=mask, affine=subject['labels'].affine),
    ) 

    return subject_new, label_names

def load_mrilab_files(path2skeleton):
    subject = tio.Subject(
            labels=tio.LabelMap(os.path.join(path2skeleton, 'labels.nii.gz')),
            T1=tio.ScalarImage(os.path.join(path2skeleton, 'T1_map.nii.gz')),
            T2s=tio.ScalarImage(os.path.join(path2skeleton, 'T2s_map.nii.gz')),
            T2=tio.ScalarImage(os.path.join(path2skeleton, 'T2_map.nii.gz')),
            PD=tio.ScalarImage(os.path.join(path2skeleton, 'Rho_map.nii.gz')),
        )
    
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

    label_names = ['Background', 'WM', 'GM', 'CSF', 'Fat', 'Skull']

    return subject, label_names

def load_custom_files(path2skeleton):
    subject = tio.Subject(
        labels=tio.LabelMap(os.path.join(path2skeleton, 'segmentation_labels.nii.gz')),
        T1=tio.ScalarImage(os.path.join(path2skeleton, 'T1_brain.nii.gz')),
        B0=tio.ScalarImage(os.path.join(path2skeleton, 'B0_brain.nii.gz')),
    )

    labels = subject['labels'].data[0]
    old_labels = labels.clone()
    # Original labelling: GM=1, WM=2, CSF=3, Bone=4, Soft=5
    # Processing to background (0), white matter (1), gray matter (2), CSF (3)
    # Fat and skull cannot be used properly because of the brain stripping for anonymization
    labels[old_labels == 1] = 2
    labels[old_labels == 2] = 1
    labels[old_labels == 3] = 3
    labels[old_labels >= 4] = 0

    # Crearte whole brain mask
    mask = (labels > 0).float()
    mask = mask.unsqueeze(0)  # Add a channel dimension
    subject['brain_mask'] = tio.LabelMap(tensor=mask, affine=subject['labels'].affine)
    
    subject['labels'].data[0] = labels

    label_names = ['Background', 'WM', 'GM', 'CSF']

    subject['labels'].affine = subject['T1'].affine  # force same affine/spacing

    return subject, label_names


def list_spectra_in_folder(folder_path, possible_components=['metabs', 'mms', 'lipids', 'water', 'noise']):
    """
    List all spectra files in a given folder.
    
    Args:
        folder_path (str): Path to the folder to search for spectra files.
    
    Returns:
        list: List of paths to spectra files.
    """
    # Holds all spectra info per folder
    spectra_by_folder = {}

    # Walk through all folders
    for root, dirs, files in os.walk(folder_path):
        nii_files = [f for f in files if f.endswith('.nii.gz')]
        json_files = [f for f in files if f.endswith('.json')]

        if nii_files:
            folder_id = os.path.relpath(root, folder_path)
            spectra_by_folder[folder_id] = {
                'components': {},
                'json': None
            }

            # Add nii.gz component files
            for file in nii_files:
                for component in possible_components:
                    if component in file:
                        spectra_by_folder[folder_id]['components'][component] = os.path.join(root, file)
                        break

            # Assume one json file per folder
            if json_files:
                spectra_by_folder[folder_id]['json'] = os.path.join(root, json_files[0])

    return spectra_by_folder


def list_invivo_spectra_in_folder(folder_path):
    """
    List all spectra files in a given folder.
    
    Args:
        folder_path (str): Path to the folder to search for spectra files.
    
    Returns:
        list: List of paths to spectra files.
    """
    # Holds all spectra info per folder
    nii_files = []
    ref_files = []

    # Walk through all folders
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith(".nii.gz"):
                must_exclude = {"H2O", "ref", "water"}
                if not any([excl in file for excl in must_exclude]):
                    nii_files.append(os.path.join(root, file))
                else:
                    ref_files.append(os.path.join(root, file))
    
    nii_files.sort()
    ref_files.sort()
    all_files = list(zip(nii_files, ref_files))
    return all_files
        