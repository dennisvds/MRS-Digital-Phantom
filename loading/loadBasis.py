####################################################################################################
#                                           loadBasis.py                                           #
####################################################################################################
#                                                                                                  #
# Authors: D.M.J. van de Sande (d.m.j.v.d.sande@tue.nl)
#          J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 30/06/22                                                                                #
#                                                                                                  #
# Purpose: Load various basis set format types for MRS and MRSI data.                              #
#                                                                                                  #
####################################################################################################


#*************#
#   imports   #
#*************#
import numpy as np
import os

from fsl_mrs.core.basis import Basis
from fsl_mrs.utils import mrs_io
from fsl_mrs.utils.misc import rescale_FID, SpecToFID, ppm2hz

from pathlib import Path

import scipy.io as sio

# own
from loading.lcmodel import read_LCModel_raw
from utils.definitions import GYRO_MAG_RATIO


#***********************************#
#   loading basis as FSL-MRS sets   #
#***********************************#
def loadBasisAsFSL(path2basis, fmt=None):
    """

    """
    if fmt is None:
        fmt = path2basis.split('.')[-1].lower()
    
    if fmt == 'basis':
        FSLbasis = mrs_io.read_basis(path2basis)
    elif fmt == 'osprey':
        FSLbasis = read_osprey_basis(path2basis)

    else:
        raise ValueError('-------------- Loading failed! --------------\n'
              ' Invalid data type or not implemented yet... \n'
              'File: ' + path2basis)
    
    return FSLbasis      

def read_osprey_basis(path2basis):
    # Load basis set
    mat_struct = sio.loadmat(path2basis)['BASIS'][0,0]

    # Get names
    name_list = list(mat_struct['name'][0])
    names = [name[0] for name in name_list]

    # Get fid array
    fid_array = mat_struct['fids'].astype(np.complex128)

    # Only keep the metabolites
    nMet = mat_struct['nMets'][0,0]
    fid_array = fid_array[:, :nMet]
    names = names[:nMet]

    # Get header information
    cf = mat_struct['Bo'][0,0].astype(float) * GYRO_MAG_RATIO['1H']
    bw = mat_struct['spectralwidth'][0,0].astype(float)
    echotime = mat_struct['te'][0,0].astype(float)
    dwelltime = mat_struct['dwelltime'][0,0].astype(float)
    fwhm = mat_struct['linewidth'][0,0].astype(float)
    headers = [{'centralFrequency': cf, 'bandwidth': bw, 'echotime': echotime, 'dwelltime': dwelltime, 'fwhm': fwhm} for _ in range(len(names))]
    
    # Apply frequency shift to fid_array
    shift_amount = 4.65 - mat_struct['centerFreq'][0,0]
    shift_amount = ppm2hz(shift_amount, cf, shift=False)
    time_axis = np.arange(fid_array.shape[0]) * dwelltime
    fid_array = fid_array * np.exp(-1j * 2 * np.pi * shift_amount * time_axis[..., np.newaxis])

    # Create basis object
    basis = Basis(fid_array, names, headers)
    
    ### Testing ###
    fids = basis.original_basis_array
    specs = np.fft.fftshift(np.fft.fft(fids, axis=0), axes=0)
    ppm = basis.original_ppm_axis + mat_struct['centerFreq'][0,0] # Osprey ppm axis
    ppm_basis = basis.original_ppm_shift_axis

    return basis