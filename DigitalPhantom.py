####################################################################################################
#                                       DigitalPhantom.py                                          #
####################################################################################################
#                                                                                                  #
# Authors: D.M.J. van de Sande (d.m.j.v.d.sande@tue.nl)                                            #
#                                                                                                  #
# Created: 26/01/2024                                                                              #
#                                                                                                  #
# Purpose: Main script to define a Digital Phantom class for MRS and MRSI simulation purposes.     #
#                                                                                                  #
#                                                                                                  #
####################################################################################################

####################################################################################################
#                                       Import packages                                            #
####################################################################################################

# Import the required packages
import os
import torch
import torchio as tio
import numpy as np
from tqdm import tqdm
import json
import nibabel as nib
from scipy.ndimage import gaussian_filter
import pandas as pd

# Import the required functions from the utils.py file
from loading.load_skeleton import load_skeleton_as_subject
from preprocessing.preprocess_phantom import create_phantom_data

####################################################################################################
#                                       DigitalPhantom                                             #
####################################################################################################
#                                                                                                  #
# Main class to define a Digital Phantom for MRS and MRSI simulation purposes.                     #
#                                                                                                  #
####################################################################################################

class DigitalPhantom():
    """
    Main class to define a Digital Phantom for MRS and MRSI simulation purposes.

    Args:
        skeleton (str): The name of the skeleton to use for the phantom (currently 'BigBrainMR' and 'MRiLab' are supported).
        phantom_resolution (int): The resolution of the phantom in mm.
        path2skeleton (str): The path to the skeleton data.
        path2phantom (str): The path to save and load the phantom data.
        path2metabs (str): The path to the metabolite data.
        from_scratch (bool): Whether to create the phantom from scratch.
        concs_std (float): The standard deviation of the metabolite concentrations. If None, the standard deviation from the metabolite dataframe is used.
        sigma_lipid_spread (float): The amount of lipid bleed in the lipid map.
        grad_metabs (list): A list of metabolites to apply a gradient to.
        grad_settings (list): A list of gradient settings [min_value, max_value, direction]. Direction can be '+x', '-x', '+y', '-y', '+z', or '-z'.

    """

    def __init__(self, skeleton='BigBrainMR', resolution=3.0, path2skeleton=None, path2phantom=None, path2metabs=None, from_scratch=False,
                 concs_std=None,
                 sigma_lipid_spread=2.0, 
                 grad_metabs=['NAA'], grad_settings=[1, 2, '+y']):
        
        print(f"Initializing Digital Phantom...")

        # Set the paths to the skeleton and metabolite data if not provided
        if path2skeleton is None:
            path2skeleton = os.path.join('data', skeleton)
        if path2metabs is None:
            path2metabs = os.path.join('data', 'metabolites')

        # Check if provided paths are valid
        if os.path.exists(path2skeleton):
            self.path2skeleton = path2skeleton
        else:
            raise ValueError(f"'{path2skeleton}' does not exist. Please provide a valid path to the skeleton data.")
        if os.path.exists(path2metabs):
            self.path2metabs = path2metabs
        else:
            raise ValueError(f"'{path2metabs}' does not exist. Please provide a valid path to the metabolite data.")
        
        # Load the metabolite dataframe
        if 'metab_df.csv' not in os.listdir(path2metabs):
            raise ValueError("metab_df.csv not found in path2metabs. Please create a metabolite dataframe first.")
        else:
            self.metab_df = pd.read_csv(os.path.join(path2metabs, 'metab_df.csv'))
        
        # Initialize the class variables
        self.skeleton = skeleton
        self.resolution = resolution
        self.from_scratch = from_scratch
        self.concs_std = concs_std
        self.sigma_lipid_spread = sigma_lipid_spread
        self.grad_metabs = grad_metabs
        self.grad_dir = grad_settings[2]
        self.grad_max = grad_settings[1]
        self.grad_min = grad_settings[0]
        # Redefine the path to the phantom based on the skeleton and resolution
        if path2phantom is None:
            raise ValueError("Please provide a path to save and/or load the phantom data.")
        self.path2phantom = os.path.join(path2phantom, f'{skeleton}', f'{resolution}mm')

        # Check if the path to the phantom exists and if directory is not empty
        # Also check if the phantom should be created from scratch
        if os.path.exists(self.path2phantom) and os.listdir(self.path2phantom) and not self.from_scratch:
            print(f'Phantom found at "{self.path2phantom}". Loading phantom...')
            self._load_phantom()
            print(f'Phantom loaded!')
        else:
            print(f'Creating a new phantom from scratch...')
            os.makedirs(self.path2phantom, exist_ok=True)
            self._create_phantom()
            self._load_phantom()
            print(f'Phantom created and saved at "{self.path2phantom}".')
    
    def get_phantom_info(self):
        print(f"###############################")
        print(f"Digital Phantom Information:")
        print(f"  - Skeleton: {self.skeleton}")
        print(f"  - Resolution: {self.resolution} mm")
        print(f"  - Metabolite Data Shape: {self.metab_data.shape}")
        print(f"  - Metabolite Data Info: {self.dim_info}")
        print(f"  - Metabolite Dimension Info: {self.metab_dim_info}")
        print(f"  - Metabolite Mapping: {self.metab_mapping}")
        print(f"  - Lipid Map Available: {self.lipid_map is not None}")
        if self.lipid_map is not None:
            print(f"    - Lipid Bleed: {self.sigma_lipid_spread}")
        print(f"  - Gradient Applied: {len(self.grad_metabs) != 0}")
        if len(self.grad_metabs) != 0:
            print(f"    - Gradient Direction: {self.grad_dir}")
            print(f"    - Gradient Max Value: {self.grad_max}")
            print(f"    - Gradient Min Value: {self.grad_min}")
        print(f"###############################")
    
    def extract_sim_data(self, metabs=[]):
        # Sort input metabolites alphabetically
        metabs = sorted(metabs, key=str.casefold)

        # Filter metabolite data based on input metabolites
        if len(metabs) != 0:
            metab_indices = [self.metab_mapping[metab] for metab in metabs]
            metab_data = self.metab_data[...,metab_indices,:]
            metab_names = metabs
        else:
            metab_data = self.metab_data
            metab_names = list(self.metab_mapping.keys())

        return metab_data, metab_names
    
    def create_metab_map(self, metab, new_resolution=None):
        metab_id = self.metab_mapping[metab]
        metab_map = self.metab_data[..., metab_id, self._find_dim_key('conc')]
        # Set nan values to zero
        metab_map[np.isnan(metab_map)] = 0
        return metab_map
    
    def _create_phantom(self):
        # Load skeleton data based on the requested skeleton
        subject = load_skeleton_as_subject(self.path2skeleton, self.skeleton, self.resolution)
        # Create phantom data
        create_phantom_data(subject, self.metab_df, self.path2phantom, save=True)
    
    def _load_phantom(self):
        # List all nii-files in path2phantom
        nii_files = [f for f in os.listdir(self.path2phantom) if f.endswith('.nii.gz')]
        subject_files = {}

        # Iterate over the files and load them
        for f in nii_files:
            # Extract the name after the first underscore and before the extension
            name = f.split('_')[1].split('.')[0]
            
            # Load the NIfTI file using nibabel
            file_path = os.path.join(self.path2phantom, f)
            nifti_data = nib.load(file_path)
            data = nifti_data.get_fdata()

            if name == 'metab':
                # Extract dimension metadata
                meta_data = nifti_data.header.extensions[0].get_content()
                metab_header = json.loads(meta_data)
                self.dim_info = metab_header['dim_info']
                self.metab_mapping = metab_header['metab_mapping']
                self.metab_dim_info = metab_header['metab_dim_info']
                self.affine = nifti_data.affine
                self.metab_data = data
            
            elif name == 'labels':
                setattr(self, name, data)
                subject_files[name] = tio.LabelMap(tensor=torch.Tensor(data).unsqueeze(0), affine=nifti_data.affine)
            
            else:
                setattr(self, name, data)
                subject_files[name] = tio.ScalarImage(tensor=torch.Tensor(data).unsqueeze(0), affine=nifti_data.affine)  

        # Create a subject object
        self.subject = tio.Subject(subject_files) 

        # Process labels, concentrations, and extract lipid map
        self.lipid_map, self.labels, self.metab_data = self._process_maps(self.labels, self.metab_data)
    
    def _process_maps(self, labels, metab_data):
        '''
        Process the labels, metabolite data, and create a lipid map if necessary.

        Args:
            labels (numpy.ndarray): A 3D numpy array representing the labels.
            metab_data (numpy.ndarray): A 4D numpy array representing the metabolite data.

        Returns:
            lipid_map (numpy.ndarray): A 3D numpy array representing the lipid map.
            labels (numpy.ndarray): A 3D numpy array representing the labels.
            metab_data (numpy.ndarray): A 4D numpy array representing the metabolite data.
        '''
        # Create lipid map
        if self.skeleton == 'BigBrainMR':
            lipid_map = None # No lipid map for BigBrain, since no information about lipid/skull region is available
        elif self.skeleton == 'MRiLab':
            lipid_map = self._generate_lipid_map(labels, lipid_region=4)
        
        # Use only background (0), white matter (1), gray matter (2), CSF (3) for labels
        labels[labels > 3] = 0

        # Apply concentration gradients if metabolites are specified
        if len(self.grad_metabs) != 0:
            gradient_map = self._generate_gradient_map(direction=self.grad_dir, max_value=self.grad_max, min_value=self.grad_min)
            metab_names = self.grad_metabs
            metab_indices = sorted([self.metab_mapping[metab] for metab in metab_names])
            gradient_map = np.repeat(gradient_map[:, :, :, np.newaxis], len(metab_indices), axis=-1)
            metab_data[...,metab_indices, self._find_dim_key('conc_mean')] = metab_data[...,metab_indices, self._find_dim_key('conc_mean')] * gradient_map

        # Set concentrations based on the mean and standard deviation
        metab_data = self._set_concentrations(std=self.concs_std)

        return lipid_map, labels, metab_data

    def _generate_gradient_map(self, direction='+y', max_value=2, min_value=1):
        """
        Generate a gradient map based on the specified direction and value range.

        Args:
            direction (str): The direction of the gradient ('+x', '-x', '+y', '-y', '+z', or '-z').
            max_value (float): The maximum value of the gradient.
            min_value (float): The minimum value of the gradient.

        Returns:
            gradient_map (numpy.ndarray): A 3D numpy array representing the gradient map.
        """
        shape = self.labels.shape
        gradient_map = np.zeros(shape)
        
        # Parse direction and sign
        direction_sign = direction[0]
        direction = direction[1:]
        
        if direction == 'x':
            for i in range(shape[0]):
                if direction_sign == '+':
                    gradient_map[i, :, :] = max(min_value, i / (shape[0] - 1) * max_value)
                elif direction_sign == '-':
                    gradient_map[i, :, :] = max(min_value, (shape[0] - 1 - i) / (shape[0] - 1) * max_value)
        elif direction == 'y':
            for j in range(shape[1]):
                if direction_sign == '+':
                    gradient_map[:, j, :] = max(min_value, j / (shape[1] - 1) * max_value)
                elif direction_sign == '-':
                    gradient_map[:, j, :] = max(min_value, (shape[1] - 1 - j) / (shape[1] - 1) * max_value)
        elif direction == 'z':
            for k in range(shape[2]):
                if direction_sign == '+':
                    gradient_map[:, :, k] = max(min_value, k / (shape[2] - 1) * max_value)
                elif direction_sign == '-':
                    gradient_map[:, :, k] = max(min_value, (shape[2] - 1 - k) / (shape[2] - 1) * max_value)
        else:
            raise ValueError("Invalid direction. Please specify '+x', '-x', '+y', '-y', '+z', or '-z'.")
        
        return gradient_map

    def _generate_lipid_map(self, labels, lipid_region=4):
        """
        Generate a lipid map using a point-spread-function.

        Args:
            labels (numpy.ndarray): A 3D numpy array representing the labels.
            lipid_region (int): The label of the lipid region.

        Returns:
            lipid_map (numpy.ndarray): A 3D numpy array representing the lipid map.
        """
        lipid_labels = (labels == lipid_region).astype(float)
        lipid_map = gaussian_filter(lipid_labels, sigma=self.sigma_lipid_spread)
        lipid_map = lipid_map / np.max(lipid_map)

        return lipid_map

    def _find_dim_key(self, key):
        for k, v in self.metab_dim_info.items():
            if v == key:
                return int(k)
        return None

    def _set_concentrations(self, std=None):
        # Find keys for conc_mean and conc_std in metab_dim_info
        concs_mean_key = self._find_dim_key('conc_mean')
        concs_std_key = self._find_dim_key('conc_std')

        if concs_mean_key is None or concs_std_key is None:
            raise ValueError("Concentration mean and standard deviation not found in the metabolite dimension info.")
        
        # Extract concs_mean and concs_std from metab_data
        concs_mean = self.metab_data[..., concs_mean_key]
        concs_std = self.metab_data[..., concs_std_key]
        
        # Sample concentrations from normal distribution
        if std is not None:
            # If external std is provided, use it for sampling
            concs = np.random.normal(concs_mean, std*concs_mean)
        else:
            # Use the internal concs_std for sampling
            concs = np.random.normal(concs_mean, concs_std)

        # Replace concs_mean with the sampled concentrations in metab_data
        self.metab_data[..., concs_mean_key] = concs
        
        # Remove concs_std from metab_data by deleting its index
        self.metab_data = np.delete(self.metab_data, concs_std_key, axis=-1)

        # Update metab_dim_info:
        # Remove concs_std, rename concs_mean to conc, and adjust the keys to reflect the new array dimensions
        new_metab_dim_info = {}
        
        shift = 0  # To adjust for the removal of the concs_std dimension
        for key, value in sorted(self.metab_dim_info.items(), key=lambda x: int(x[0])):
            key = int(key)  # Convert string key to integer for logic comparison
            if value == 'conc_std':
                shift = -1  # Once conc_std is removed, shift the subsequent indices
            elif value == 'conc_mean':
                new_metab_dim_info[str(key)] = 'conc'  # Rename conc_mean to conc
            else:
                new_metab_dim_info[str(key + shift)] = value  # Adjust keys for other entries

        # Assign the updated dim info back
        self.metab_dim_info = new_metab_dim_info

        return self.metab_data



        
        




