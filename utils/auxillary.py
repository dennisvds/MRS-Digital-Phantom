import re
import numpy as np
import nibabel as nib
import json
from scipy.fft import fft, fftshift, ifft, ifftshift
from tqdm import tqdm
from scipy.signal import hilbert
import time
from contextlib import contextmanager
import os
import h5py
from scipy.ndimage import zoom
import dask.array as da
from dask.diagnostics import ProgressBar

def extract_todos(filepath):
    todos = []
    with open(filepath, 'r') as file:
        lines = file.readlines()
        for line in lines:
            todo_match = re.search(r'#\s*TODO:\s*(.*)', line)
            if todo_match:
                todos.append(todo_match.group(1))
    return todos

@contextmanager
def timer():
    start = time.time()
    yield
    end = time.time()
    print(f"Duration: {end - start:.4f} seconds")

def lorentzian(ppm, amplitude, fwhm, center, phase, offset):
    spec = np.sqrt(2/np.pi) * (fwhm/2 - 1j * (ppm - center)) / (fwhm**2 + (ppm-center)**2) * amplitude * np.exp(1j * phase * np.pi/180)
    spec = spec + offset
    return spec.real


def save_nifti_mrsi(spec_data, basis, affine, etime, rtime, path2save, save_name):
    # Transfrom spectral data to time domain
    # spec_data shape: (x, y, z, spectral_points)
    print("Saving mrs-nifti file...")
    fid_data = spec2fid(spec_data, axis=-1)
    fid_data = np.conjugate(fid_data)

    metadata = {'SpectrometerFrequency': [basis.cf],
                'ResonantNucleus': ['1H',], 
                'RepetitionTime': rtime,
                'EchoTime': etime*1e-3, 
                'Manufacturer': 'Synthetic' }    
    
    # dimensions_dict = {'dim_5': 'DIM_DYN'}   
    
    full_metadata = metadata
    # full_metadata.update(dimensions_dict)
    metadata_json = json.dumps(full_metadata)

    # Create nifti object
    newobj = nib.nifti2.Nifti2Image(fid_data, affine=affine)

    # Write new header
    pixDim = newobj.header['pixdim']
    pixDim[1:4] = affine.diagonal()[:3]
    pixDim[4] = basis.dwelltime
    newobj.header['pixdim'] = pixDim

    # Set q_form = 0
    newobj.header.set_qform(None,code=0)

    # Set conformance level 
    newobj.header['intent_name'] = b'mrs_v0_2'
    newobj.header.extensions.clear()

    # Write extension
    extension = nib.nifti1.Nifti1Extension(44, metadata_json.encode('UTF-8'))
    newobj.header.extensions.append(extension)

    # # From nii obj and write 
    resolution = np.diag(affine)[:3]
    nib.save(newobj, path2save + f'/{save_name}_{resolution[0]}_{resolution[1]}_{resolution[2]}mm.nii.gz')

    print(f"Done!")

def save_numpy_mrsi(mrsi_data, save_dir, save_name, phantom):
    # Check if save_dir exists
    if not os.path.exists(save_dir):
        raise FileNotFoundError(f"Directory '{save_dir}' does not exist.")
    # Create save path
    save_path = os.path.join(save_dir, phantom.skeleton, str(phantom.resolution)+'mm', save_name, '.npy')
    # Create directory if it does not exist
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    # Save numpy file
    print(f"Saving numpy file to {save_path}...")
    np.save(save_path, mrsi_data)
    print(f"Saving numpy file to {save_path} done!")

    return save_path

def save_hdf5_mrsi(mrsi_data, save_dir, save_name, phantom):
    # Check if save_dir exists
    if not os.path.exists(save_dir):
        raise FileNotFoundError(f"Directory '{save_dir}' does not exist.")
    # Create save path
    save_path = os.path.join(save_dir, phantom.skeleton, f"{phantom.resolution}mm", f"{save_name}.h5")
    # Create directory if it does not exist
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    
    # Save HDF5 file
    print(f"Saving MRSI data as HDF5 file to {save_path}...")
    with h5py.File(save_path, 'w') as h5f:
        h5f.create_dataset('mrsi_data', data=mrsi_data)
    print(f"Saving HDF5 file done!")
    
    return save_path


def fid2spec(fid_data, axis=-1):
    return fftshift(fft(fid_data, axis=axis), axes=axis)

def spec2fid(spec_data, axis=-1):
    return ifft(ifftshift(spec_data, axes=axis), axis=axis)

def downsample_mrsi(filename, original_affine, target_resolution):
    print('Starting downsampling...')
    with h5py.File(filename, 'r') as h5f:
        data = da.from_array(h5f['mrsi_data'], chunks=(32, 32, 32, -1))
        original_resolution = np.diag(original_affine)[:3]

        # Calculate block_size
        target_resolution = np.array(target_resolution)
        block_size = np.ceil(target_resolution / original_resolution).astype(int)
        new_resolution = original_resolution * block_size

        # Check if target resolution is lower than original resolution
        if not np.all(target_resolution >= original_resolution):
            raise ValueError(f"Target resolution must be lower than original resolution. Original resolution: {tuple(original_resolution)} mm, target resolution: {tuple(target_resolution)} mm")

        # Calculate the required padding
        pad_width = [
            (0, (block_size[i] - data.shape[i] % block_size[i]) % block_size[i])
            for i in range(3)
        ] + [(0, 0)]  # No padding for spectral axis

        # Notify if padding is needed
        if any(pad[0] + pad[1] > 0 for pad in pad_width):  # Only print if padding is non-zero
            print(f"Padding MRSI data with {pad_width}...")

        # Pad the data to make dimensions divisible by block_size
        padded_data = da.pad(data, pad_width, mode='constant', constant_values=0.0)

        # Downsample using block averaging
        block_shape = {i: block_size[i] for i in range(3)}  # Spatial axes (0, 1, 2)
        block_shape[3] = 1  # Spectral axis (unchanged)
        reduced_data = da.coarsen(np.mean, padded_data, block_shape)

        # Compute the reduced data
        print(f"Downsampling metabolic map from {tuple(original_resolution)} mm to {tuple(new_resolution)} mm using block averaging...")
        with ProgressBar():
            reduced_data = reduced_data.compute()

        # Remove padding from the reduced data
        crop_slices = tuple(
            slice(0, (data.shape[i] // block_size[i]) * block_size[i] // block_size[i])
            for i in range(3)
        ) + (slice(None),)  # Keep spectral axis unchanged
        reduced_data = reduced_data[crop_slices]

    new_affine = np.copy(original_affine)
    new_affine[:3, :3] = np.diag(new_resolution)
    print(f"Downsampling completed!")
    return reduced_data, new_affine

def downsample_metab_map(metab_map, original_affine, target_resolution):
    print("Starting metabolic map downsampling...")
    # Convert NumPy array to Dask array, with adjusted chunk size
    metab_map = da.from_array(metab_map, chunks=(32, 32, 32))  # You can tweak the chunk size

    original_resolution = np.diag(original_affine)[:3]

    # Calculate block_size based on target resolution
    target_resolution = np.array(target_resolution)
    block_size = np.ceil(target_resolution / original_resolution).astype(int)
    new_resolution = original_resolution * block_size

    # Check if target resolution is lower than original resolution
    if not np.all(target_resolution >= original_resolution):
        raise ValueError(f"Target resolution must be lower than original resolution. Original resolution: {tuple(original_resolution)} mm, target resolution: {tuple(target_resolution)} mm")

    # Calculate new shape after downsampling
    new_shape = np.ceil(np.array(metab_map.shape) / block_size).astype(int)

    # Calculate padding needed for each axis to ensure the dimensions are divisible by block_size
    pad_width = [(0, (block_size[i] - (metab_map.shape[i] % block_size[i])) % block_size[i]) for i in range(3)]
    
    # Notify if padding is needed
    if any(pad[0] + pad[1] > 0 for pad in pad_width):  # Only print if padding is non-zero
        print(f"Padding metabolic map with {pad_width}...")
    
    # Pad the metab_map
    padded_map = da.pad(metab_map, pad_width, mode='constant', constant_values=0.0)

    # Define block shape for coarsening (downsampling)
    block_shape = {i: block_size[i] for i in range(3)}  # Only for spatial dimensions

    # Downsample using block averaging
    print(f"Downsampling metabolic map from {tuple(original_resolution)} mm to {tuple(new_resolution)} mm using block averaging...")
    reduced_map = da.coarsen(np.mean, padded_map, block_shape)

    # Compute the reduced map (lazy evaluation until this point)
    # Add progress bar to monitor computation
    with ProgressBar():
        reduced_map = reduced_map.compute()

    # Update affine transformation for new resolution
    new_affine = np.copy(original_affine)
    new_affine[:3, :3] = np.diag(new_resolution)

    print("Metabolic map downsampling completed!")
    return reduced_map, new_affine




def get_class_info(cls):
    class_variables = []
    instance_variables = []
    methods = []
    
    for item in dir(cls):
        # Exclude built-in methods and class methods
        if item.startswith("__") and item.endswith("__"):
            continue
        
        value = getattr(cls, item)
        
        # Check if it's a class variable and not starting with '_'
        if not callable(value) and not isinstance(value, staticmethod) and not item.startswith('_'):
            class_variables.append(item)
            
        # Check if it's an instance variable and not starting with '_'
        elif not callable(value) and not isinstance(value, classmethod) and not isinstance(value, staticmethod) and not item.startswith('_'):
            instance_variables.append(item)
            
        # Check if it's a method and not starting with '_'
        elif callable(value) and not isinstance(value, classmethod) and not isinstance(value, staticmethod) and not item.startswith('_'):
            methods.append(item)
    
    return class_variables, instance_variables, methods


def scale_T1_map(t1_map, labels, tesla=3.0):
    # Scaling factors from Rooney et al. (2007), doi:10.1002/mrm.21122

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


#*********************************#
#   create random walk artifact   #
#*********************************#
def randomWalk(ppm_axis, rw_range=[4.5, 5.5], scale=1, smooth=10, ylim=[-1, 1]):
    """
        Produces a spectrum with a random walk within a specified ppm range.
        
        @param ppm_axis -- The ppm axis of the spectrum (an array of ppm values).
        @param rw_range -- The range in ppm where the random walk is applied.
        @param scale -- The y scale of the steps of the walk.
        @param smooth -- The smoothness of the walk.
        @param ylim -- The y limits, format [min, max].

        @returns -- The spectrum with the random walk applied only in the rw_range.
    """
    waveLength = len(ppm_axis)  # length of the ppm axis
    spectrum = np.zeros(waveLength)  # initialize spectrum with zeros

    # Find the indices corresponding to the rw_range within the ppm_axis
    start_idx = np.abs(ppm_axis - rw_range[0]).argmin()
    end_idx = np.abs(ppm_axis - rw_range[1]).argmin()

    # Initialize the random walk in the rw_range
    y = np.random.uniform(ylim[0], ylim[1])  # init randomly between limits
    wave = []
    for _ in range(end_idx - start_idx):
        step = np.random.normal(scale=scale)
        if ylim[0] <= y + step <= ylim[1]:
            y += step
        else:
            y -= step
        wave.append(y)
    
    # Insert the random walk into the spectrum at the correct indices
    spectrum[start_idx:end_idx] = wave

    # Smooth the spectrum
    spectrum = np.convolve(spectrum, np.ones(smooth) / smooth, mode='same')

    # Put in the complex domain
    spectrum = hilbert(spectrum)

    return spectrum

#*********************************#
#   create random peak artifact   #
#*********************************#
def randomPeak(waveLength=1024, batch=1, amp=None, pos=None, width=None, phase=None, td=False):
    """
        Produces a spectrum of a random peak.

        @param waveLength -- The number of points of the spectrum.
        @param batch -- The number of peaks to produce.
        @param amp -- The amplitude of the peak.
        @param pos -- The position of the peak.
        @param width -- The width of the peak.
        @param phase -- The phase of the peak.
        @param td -- If True the spectrum is returned in the time domain.

        @returns -- The random wave shape from the random peak (complex).
    """
    if amp is None: amp = np.ones((batch, 1))
    if pos is None: pos = np.ones((batch, 1)) * waveLength // 2
    if width is None: width = np.ones((batch, 1)) * 10
    if phase is None: phase = np.zeros((batch, 1)) + 0.5
    t = np.arange(waveLength)[None, :]

    x = amp * np.exp(- (t - pos) ** 2 / (2 * width ** 2))
    x = hilbert(x.real) * np.exp(- 1j * phase)

    if td: x = np.fft.ifft(x, axis=-1)
    return x


def baseline_init(points, order, first, last):
    '''
    Create a baseline regressor for the MRS signal.
    The regressor is a polynomial of order 'order' with complex coefficients.

    @param points -- The number of points in the spectrum.
    @param order -- The order of the polynomial.
    @param first -- The first index of the regressor.
    @param last -- The last index of the regressor.

    @returns -- The baseline regressor.
    '''
    x = np.zeros(points, complex)
    x[first:last] = np.linspace(-1, 1, last - first)
    B = []
    for i in range(order + 1):
        regressor = x ** i
        if i > 0:
            regressor = regress_out(regressor, B, keep_mean=False)

        B.append(regressor.flatten())
        B.append(1j * regressor.flatten())

    B = np.asarray(B).T
    tmp = B.copy()
    B = 0 * B
    B[first:last, :] = tmp[first:last, :].copy()
    return B

#***************#
#   regressor   #
#***************#
def regress_out(x, conf, keep_mean=True):
    """
    Linear deconfounding

    Ref: Clarke WT, Stagg CJ, Jbabdi S. FSL-MRS: An end-to-end spectroscopy analysis package.
    Magnetic Resonance in Medicine 2021;85:2950–2964 doi: https://doi.org/10.1002/mrm.28630.
    """
    if isinstance(conf, list):
        confa = np.squeeze(np.asarray(conf)).T
    else:
        confa = conf
    if keep_mean:
        m = np.mean(x, axis=0)
    else:
        m = 0
    return x - confa @ (np.linalg.pinv(confa) @ x) + m



