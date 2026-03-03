from fsl_mrs.utils.preproc import nifti_mrs_proc as proc
import fsl_mrs.utils.mrs_io as mrs_io
from nifti_mrs.create_nmrs import gen_nifti_mrs
import nibabel as nib
import json
import numpy as np
import os
import random
import datetime
import string
from nibabel.nifti1 import Nifti1Extension

from utils.auxillary import spec2fid, fid2spec


def save_nifti_mrs(save_dir, total_selection, individual_selection, components_data, affine, dwelltime,
                   spec_freq, nucleus='1H', sim_params=None, base_name=None):
    """
    Save selected MRS components as NIfTI-MRS files, each in its own subfolder.

    Parameters:
    - save_dir (str): Directory where output should be stored.
    - total_selection (list of str): Component names to combine and save as "total".
    - individual_selection (list of str): Component names to save individually.
    - components_data (dict): Mapping from component names to complex data arrays.
    - affine (np.ndarray): 4x4 affine matrix for orientation.
    - dwelltime (float): Dwell time (s).
    - spec_freq (float): Spectrometer frequency (MHz).
    - nucleus (str): Resonant nucleus, default '1H'.
    - sim_params (dict): Optional simulation parameter metadata.
    - base_name (str): Optional base name for output folder.
    """

    def _generate_base_name():
        """Create a unique name for the phantom."""
        date_str = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        rand_id = ''.join(random.choices(string.ascii_lowercase + string.digits, k=4))
        return f"phantom_{date_str}_{rand_id}"

    def _prepare_data(data):
        """Transform spectrum to FID and ensure shape is 4D."""
        data = spec2fid(data, axis=-1)
        if data.ndim == 1:
            data = data[np.newaxis, np.newaxis, np.newaxis, :, np.newaxis]  # Add dimensions for 4D
            data = np.conjugate(data)  # Conjugate to match NIfTI-MRS convention
        return np.array(data, dtype=np.complex64)

    def _create_metadata():
        """Assemble NIfTI-MRS metadata dictionary."""
        dimensions = {'dim_5': 'DIM_DYN'}
        metadata = {
            "SpectrometerFrequency": [float(spec_freq)],
            "ResonantNucleus": [nucleus],
            "RepetitionTime": float(repetition_time or 0),
            "EchoTime": float(echotime or 0),
            "Manufacturer": manufacturer,
            # "DwellTime": float(dwelltime),
            # "Global": {}
        }
        full_metadata = metadata.copy()
        full_metadata.update(dimensions)
        return full_metadata

    def _save_nifti_mrs_image(data, metadata, path):
        """Save a NIfTI-MRS image with header and extension metadata."""
        img = nib.nifti1.Nifti1Image(data, affine=affine)
        pixdim = img.header['pixdim']
        pixdim[1:4] = affine.diagonal()[:3]
        pixdim[4] = dwelltime
        img.header['pixdim'] = pixdim
        img.header.set_qform(None, code=0)
        img.header['intent_name'] = b'mrs_v0_2'
        img.header.extensions.clear()
        extension = nib.nifti1.Nifti1Extension(44, json.dumps(metadata).encode('utf-8'))
        img.header.extensions.append(extension)
        nib.save(img, path)

    # Determine output folder
    if base_name is None:
        base_name = _generate_base_name()
    spectrum_folder = os.path.join(save_dir, base_name)
    os.makedirs(spectrum_folder, exist_ok=True)

    # Extract sim metadata if present
    echotime = repetition_time = None
    manufacturer = 'Unknown'
    if sim_params:
        bs = sim_params.get('basis_set_settings', {})
        sp = sim_params.get('simulation_params', {})
        echotime = bs.get('TE', 0) * 1e-3  # ms to s
        repetition_time = sp.get('TR', 0) * 1e-3
        manufacturer = sp.get('manufacturer', 'Unknown')

    # Save total spectrum
    if total_selection:
        total_data = sum(components_data[name] for name in total_selection)
        total_data = _prepare_data(total_data)
        metadata = _create_metadata()
        total_path = os.path.join(spectrum_folder, "total.nii.gz")
        _save_nifti_mrs_image(total_data, metadata, total_path)

    # Save individual components (using nibabel same as total)
    for name in individual_selection:
        data = _prepare_data(components_data[name])
        metadata = _create_metadata()
        safe_name = name.replace(' ', '_').lower()
        comp_path = os.path.join(spectrum_folder, f"{safe_name}.nii.gz")
        _save_nifti_mrs_image(data, metadata, comp_path)

    # Save simulation parameters if available
    if sim_params:
        with open(os.path.join(spectrum_folder, "sim_params.json"), 'w') as f:
            json.dump(sim_params, f, indent=4)

    print(f"Saved spectrum to {spectrum_folder}")


def load_nifti_mrs(file_path):
    """
    Load a NIfTI-MRS file and return the data and metadata.
    
    Args:
        file_path (str): Path to the NIfTI-MRS file.
    
    Returns:
        tuple: A tuple containing the data and metadata.
    """
    img = nib.load(file_path)

    data = img.get_fdata(dtype=np.complex64)
    header = img.header
    dwelltime = header['pixdim'][4]


    hdr_ext_codes = img.header.extensions.get_codes()
    mrs_hdr_ext = json.loads(img.header.extensions[hdr_ext_codes.index(44)].get_content())
    
    spectrometer_frequency = mrs_hdr_ext['SpectrometerFrequency'][0]
    bandwidth = 1/dwelltime
    frequency_axis = np.linspace(bandwidth/2, -bandwidth/2, data.shape[3])
    ppm_axis = 4.65+frequency_axis/spectrometer_frequency

    return ppm_axis, data[0,0,0,:], mrs_hdr_ext


def preprocess_press(mrs_file_pair, basis_file=None):
    # Preprocessing according to: https://open.win.ox.ac.uk/pages/fslcourse/practicals/fsl_mrs/index.html
    # and: https://git.fmrib.ox.ac.uk/saad/fsl_mrs/-/blob/master/example_usage/Example%20SVS%20processing%20-%20interactive%20notebook.ipynb

    # report_dir = "./reports"
    report_dir = None

    raw_act = mrs_io.read_FID(mrs_file_pair[0])
    raw_ref= mrs_io.read_FID(mrs_file_pair[1])

    # Set processing data
    act_data = raw_act
    ref_data = raw_ref

    # Coil combination
    if 'DIM_DYN' in ref_data.dim_tags:
        avg_ref_data = proc.average(ref_data, 'DIM_DYN', figure=False, report=report_dir)
    else:
        avg_ref_data = ref_data
    
    if 'DIM_COIL' in act_data.dim_tags:
        act_data = proc.coilcombine(act_data, reference=avg_ref_data, noise=None, figure=False, report=report_dir)
        ref_data = proc.coilcombine(ref_data, reference=avg_ref_data, noise=None, figure=False, report=report_dir)


    # Align dynamics (1st iteration)
    if 'DIM_DYN' in act_data.dim_tags:
        act_data = proc.align(act_data, 'DIM_DYN', ppmlim=(0.2, 4.2), figure=False, report=report_dir)

    # Remove bad quality spectra
    if 'DIM_DYN' in act_data.dim_tags:
        act_data, bad_data = proc.remove_unlike(act_data, sdlimit=2.58, niter=1, ppmlim=(0, 4.2), figure=False, report=report_dir)
        if bad_data is None:
            bad_shape = 0
        elif bad_data.ndim == 4:
            bad_shape = 1
        else:
            bad_shape = bad_data.shape[4]
        print(f"{bad_shape} bad averages identified and removed...")

    # Align Dynamics (2nd iteration)
    if 'DIM_DYN' in act_data.dim_tags:
        act_data = proc.align(act_data, 'DIM_DYN', ppmlim=(0.2, 4.2), figure=False, report=report_dir)
    if 'DIM_DYN' in ref_data.dim_tags:
        ref_data = proc.align(ref_data, 'DIM_DYN', ppmlim=(0, 8), figure=False, report=report_dir)

    # Averaging
    if 'DIM_DYN' in act_data.dim_tags:
        act_data = proc.average(act_data, 'DIM_DYN', figure=False, report=report_dir)
    if 'DIM_DYN' in ref_data.dim_tags:
        ref_data = proc.average(ref_data, 'DIM_DYN', figure=False, report=report_dir)

    # # Eddy current correction
    # if 'DIM_DYN' in ref_data.dim_tags:
    #     eccRef = proc.average(ref_data, 'DIM_DYN')
    # else:
    #     eccRef = ref_data

    # act_data = proc.ecc(act_data, eccRef, figure=False, report=report_dir)
    # ref_data = proc.ecc(ref_data, eccRef, figure=False, report=report_dir)

    # Leftshift/Truncation (only if needed)
    # act_data = proc.truncate_or_pad(act_data, -1, 'first', figure=False, report=report_dir)
    # ref_data = proc.truncate_or_pad(ref_data, -1, 'first', figure=False, report=report_dir)

    # Residual water removal
    hlsvdlimits = [-0.25, 0.25]
    act_data = proc.remove_peaks(act_data, hlsvdlimits, limit_units='ppm', figure=False, report=report_dir)

    # Reference shift to Cr (3.027 ppm)
    act_data = proc.shift_to_reference(act_data, 3.027, (2.9, 3.1), figure=False, report=report_dir)

    # Apply phasing on a single peak (tCr)
    act_data = proc.phase_correct(act_data, (2.9, 3.1), figure=False, report=report_dir)

    final_wref = proc.phase_correct(ref_data, (4.55, 4.7), hlsvd=False, figure=False, report=report_dir)

    # Apply apodization
    # act_data = proc.apodize(act_data, amount=[10.0], figure=False, report=report_dir)

    # Generate MRS object
    mrs=act_data.mrs(basis_file=basis_file, ref_data=final_wref)
    
    if report_dir is not None:
        print("Create preprocessing report...")
        description = "output"
        out_folder = report_dir
        file_name = f"report_{os.path.basename(mrs_file_pair[0]).rsplit('.')[0]}.html"
        delete_name = out_folder + "/*.html"
        os.system(f"merge_mrs_reports -d {description} -o {out_folder} -f {file_name} --delete {delete_name}")


    return mrs

def preprocess_press_sim(sim_file):
    
    report_dir = None
    act_data = mrs_io.read_FID(sim_file)

    # Residual water removal
    hlsvdlimits = [-0.25, 0.25]
    act_data = proc.remove_peaks(act_data, hlsvdlimits, limit_units='ppm', figure=False, report=report_dir)

    # Reference shift to Cr (3.027 ppm)
    act_data = proc.shift_to_reference(act_data, 3.027, (2.9, 3.1), figure=False, report=report_dir)

    # Apply phasing on a single peak (tCr)
    act_data = proc.phase_correct(act_data, (2.9, 3.1), figure=False, report=report_dir)

    # Generate MRS object
    mrs=act_data.mrs(basis_file=None)

    if report_dir is not None:
        print("Create preprocessing report...")
        description = "output"
        out_folder = report_dir
        file_name = f"report_{os.path.basename(sim_file).rsplit('.')[0]}.html"
        delete_name = out_folder + "/*.html"
        os.system(f"merge_mrs_reports -d {description} -o {out_folder} -f {file_name} --delete {delete_name}")
    return mrs

