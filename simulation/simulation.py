####################################################################################################
#                                        simulation.py                                             #
####################################################################################################
#                                                                                                  #
# Authors: D.M.J. van de Sande (d.m.j.v.d.sande@tue.nl)                                            #
#                                                                                                  #
# Created: 03/02/23                                                                                #
#                                                                                                  #
# Purpose: Defines the main methods to generate batches of spectra using the Digital MRS Phantom.  #
#                                                                                                  #
####################################################################################################

#*************#
#   imports   #
#*************#
import numpy as np
from tqdm import tqdm
from scipy.signal import hilbert


def simulate_spectra(phantom, basis, sigModel, metabs=None, slice_range=[],  batch_size=500):
    # Prepare simulation data
    data, coordinates, metab_names, metab_dim_info = prepare_sim_data(phantom, basis, metabs, slice_range)

    # Initialize output array
    out = np.zeros((*phantom.metab_data.shape[0:3], basis.n), dtype=np.complex64)

    print(f"Simulating spectra for {len(data)} voxels...")
    for i in tqdm(range(0, len(data), batch_size)):
        # Get batches of data
        coordinates_batch = coordinates[i:i + batch_size]
        data_batch = data[i:i + batch_size]
        labels_batch = phantom.labels[coordinates_batch[:, 0], coordinates_batch[:, 1], coordinates_batch[:, 2]]
        try:
            lipid_mask_batch = phantom.lipid_map[coordinates_batch[:, 0], coordinates_batch[:, 1], coordinates_batch[:, 2]]
        except:
            lipid_mask_batch = None

        # Simulate spectra
        spectra = sigModel.forward(data_batch, labels_batch, lipid_mask_batch)

        # Store simulated spectra
        for j, (x, y, z) in enumerate(coordinates_batch):
            out[x, y, z] = spectra[j]   
    
    # Get ppm axis 
    ppm_axis = basis.ppm
    print("Spectral simulation completed!")
    
    return ppm_axis, out


def prepare_sim_data(phantom, basis, metabs=None, slice_range=[]):
    data, metab_names = phantom.extract_sim_data(metabs)
    metab_dim_info = phantom.metab_dim_info
    labels = phantom.labels

    # Filter the data if certain metabolites are not present in the basis
    data, metab_names = filter_metabolite_names(data, metab_names, basis)

    # Get data shape (x, y, and z dimensions)
    data_shape = data.shape[0:3]

    # Get x, y, z coordinates
    x, y, z = np.meshgrid(np.arange(data_shape[0]), np.arange(data_shape[1]), np.arange(data_shape[2]), indexing='ij')

    # Select all slices if slice_range is not provided
    if len(slice_range) == 0:
        slice_range = [0, data_shape[2]]

    # Prepare masks for voxels that need to be simulated
    z_slices = np.arange(max(0, slice_range[0]), min(slice_range[1], data.shape[2]))
    z_mask = np.isin(z, z_slices)

    # Get non-background voxels and select specific slice(s)
    x_nb, y_nb, z_nb = np.where(labels != 0)                        # Filter out background voxels
    x_nb, y_nb, z_nb = np.where((labels != 0) & z_mask)             # Select specific slice(s)

    # Filter out background spectra and select specific slice(s)
    data = data[x_nb, y_nb, z_nb]
    coordinates = np.stack([x_nb, y_nb, z_nb], axis=-1)

    return data, coordinates, metab_names, metab_dim_info

def filter_metabolite_names(data, metabolite_names, basis):
    # Get metabolite names from basis
    basis_metab_names = basis.names

    # Get metabolites that are not present in both the basis and the data
    metabs_to_remove = [m for m in metabolite_names if m not in basis_metab_names]

    # Remove metabolites that are not present in the basis
    if len(metabs_to_remove) > 0:
        print(f"There are metabolites in the data that are not present in the basis: {metabs_to_remove}")
        print(f"Removing these metabolites from the data...")
        metab_indices = [metabolite_names.index(m) for m in metabs_to_remove]
        data = np.delete(data, metab_indices, axis=-2)
        metabolite_names = [m for m in metabolite_names if m not in metabs_to_remove]

    return data, metabolite_names