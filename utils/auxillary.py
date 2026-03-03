import os
import pandas as pd
import torch
from scipy.fft import fft, fftshift, ifft, ifftshift
from torch.fft import fft as torch_fft, ifft as torch_ifft, fftshift as torch_fftshift, ifftshift as torch_ifftshift
import numpy as np
import nibabel as nib
import datetime
from nifti_mrs.create_nmrs import gen_nifti_mrs
import json
import random
import string
import itertools
from collections import defaultdict
import ast

from fsl_mrs.utils.preproc import nifti_mrs_proc as proc
import fsl_mrs.utils.mrs_io as mrs_io


def fid2spec(fid_data, axis=-1, shift=True):
    if isinstance(fid_data, np.ndarray):
        if shift:
            return fftshift(fft(fid_data, axis=axis), axes=axis)
        else:
            return fft(fid_data, axis=axis)
    elif isinstance(fid_data, torch.Tensor):
        if shift:
            return torch_fftshift(torch_fft(fid_data, dim=axis), dim=axis)
        else:
            return torch_fft(fid_data, dim=axis)
    else:
        raise ValueError("Input data type not supported. Must be either numpy array or torch tensor.")

def spec2fid(spec_data, axis=-1, shift=True):
    if isinstance(spec_data, np.ndarray):
        if shift:
            return ifft(ifftshift(spec_data, axes=axis), axis=axis)
        else:
            return ifft(spec_data, axis=axis)
        
    elif isinstance(spec_data, torch.Tensor):
        if shift:
            return torch_ifft(torch_ifftshift(spec_data, dim=axis), dim=axis)
        else:
            return torch_ifft(spec_data, dim=axis)
    else:
        raise ValueError("Input data type not supported. Must be either numpy array or torch tensor.")
    
def calc_scale_factor(area_1, num_protons_1, area_2, num_protons_2):
    '''
    Calculate the scaling factor between two peaks using the area under the curve and the number of protons.

    @param area_1 -- The area of the first peak.
    @param num_protons_1 -- The number of protons of the first peak.
    @param area_2 -- The area of the second peak.
    @param num_protons_2 -- The number of protons of the second peak.

    @returns -- The scaling factor. Can be used to scale the second peak to the first peak. (e.g. peak_2 * scale_factor = peak_1)
    '''

    scale_factor = (area_1 / num_protons_1) / (area_2 / num_protons_2)
    return scale_factor

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


def split_config_into_simulations(config_path):
    """
    Splits a single config file into multiple simulation_config_X.json files,
    based on the list parameters inside the config.
    
    Special case:
    - The 'metabs' parameter must be provided as a list of lists if multiple metabolite sets are intended.
    
    Args:
        config_path (str): Path to the base configuration JSON file.
    """
    # === Load config ===
    with open(config_path, 'r') as f:
        base_config = json.load(f)
    
    # === Find list parameters ===
    list_params = {}

    def find_list_params(config, parent_key=""):
        """Recursively finds parameters that are lists."""
        for key, value in config.items():
            if isinstance(value, list):
                list_params[parent_key + key] = value
            elif isinstance(value, dict):
                find_list_params(value, parent_key + key + ".")

    find_list_params(base_config)
    
    # === Prepare combinations ===
    if isinstance(base_config.get("metabs", None), list):
        metab_combinations = base_config["metabs"]
    else:
        metab_combinations = []
    
    param_combinations = []
    for param, values in list_params.items():
        if param == "metabs":
            param_combinations.append(metab_combinations)
        else:
            param_combinations.append(values)

    all_param_combinations = list(itertools.product(*param_combinations))

    # === Save each configuration ===
    output_dir = os.path.dirname(config_path)
    os.makedirs(output_dir, exist_ok=True)

    for index, param_comb in enumerate(all_param_combinations):
        sim_config = json.loads(json.dumps(base_config))  # Deep copy
        
        param_index = 0
        for param in list_params:
            value = param_comb[param_index]
            keys = param.split(".")
            temp_config = sim_config
            for key in keys[:-1]:
                temp_config = temp_config.get(key, {})
            temp_config[keys[-1]] = value
            param_index += 1

        output_file = os.path.join(output_dir, f"simulation_config_{index}.json")
        with open(output_file, 'w') as f:
            json.dump(sim_config, f, indent=4)
        print(f"Saved config {output_file}")

    print(f"Generated {len(all_param_combinations)} simulation configuration files.")

def group_simulation_configs_by_settings(config_path):
    config_dir = os.path.dirname(config_path)
    all_files = [
        os.path.join(config_dir, f)
        for f in os.listdir(config_dir)
        if 'simulation' in f and f.endswith('.json')
    ]

    groups = defaultdict(list)

    def extract_settings_key(file_path):
        with open(file_path, 'r') as f:
            config = json.load(f)
        
        # Remove voxel_definitions for grouping
        config_copy = dict(config)
        config_copy.pop("voxel_definitions", None)

        # Make a sorted tuple of (key, value) pairs for hashing
        return tuple(sorted(flatten_dict(config_copy).items()))
    
    def flatten_dict(d, parent_key=''):
        """Flattens a nested dictionary into a single level."""
        items = []
        for k, v in d.items():
            new_key = parent_key + '.' + k if parent_key else k
            if isinstance(v, dict):
                items.extend(flatten_dict(v, new_key).items())
            else:
                items.append((new_key, str(v)))
        return dict(items)

    for file in all_files:
        key = extract_settings_key(file)
        groups[key].append(file)
    
    return groups

def try_parse_value(value):
    """Tries to parse a value, converting strings like '30' to integers, '[1, 2, 3]' to lists, etc."""
    try:
        # If the value is a string that looks like a list or number, parse it
        parsed_value = ast.literal_eval(value)
        return parsed_value
    except (ValueError, SyntaxError):
        # If it's not parsable, return it as is (e.g., it's a string or other non-parsable value)
        return value

def unflatten_dict(settings_key):
    """Takes a tuple of key-value pairs (from settings_key) and rebuilds a nested dictionary with automatic parsing."""
    nested = {}

    # Step 1: Convert settings_key (tuple) into a flat dict
    if isinstance(settings_key, tuple):
        flat_dict = dict(settings_key)
    else:
        flat_dict = settings_key

    # Step 2: Unflatten
    for compound_key, value in flat_dict.items():
        keys = compound_key.split('.')
        d = nested
        for key in keys[:-1]:
            if key.isdigit():
                key = int(key)  # Interpret numbers as list indices
            if key not in d:
                d[key] = {} if not isinstance(key, int) else []
            d = d[key]
        
        # Last key handling and parsing the value
        last_key = keys[-1]
        if last_key.isdigit():
            last_key = int(last_key)

        parsed_value = try_parse_value(value)
        d[last_key] = parsed_value

    return nested

