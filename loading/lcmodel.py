####################################################################################################
#                                            lcmodel.py                                            #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#                                                                                                  #
# Created: 30/06/22                                                                                #
#                                                                                                  #
# Purpose: Load the raw and h2o file fomats of LCModel.                                            #
#                                                                                                  #
####################################################################################################


#*************#
#   imports   #
#*************#
import numpy as np


#*************************#
#   loading LCModel raw   #
#*************************#
def read_LCModel_raw(filename, conjugate=True):
    """
    Read LCModel (.RAW, .raw, and .H2O) file format. Adapted from [1].

    [1] Clarke, W.T., Stagg, C.J., and Jbabdi, S. (2020). FSL-MRS: An end-to-end
        spectroscopy analysis package. Magnetic Resonance in Medicine, 85, 2950 - 2964.

    @param filename -- Path to .RAW/.H2O file.
    @param bool conjugate -- Apply conjugation upon read.

    @returns -- The basis set data/FID and header if possible.
    """
    header = []
    data   = []
    in_header = False
    after_header = False
    with open(filename, 'r') as f:
        for line in f:
            if (line.find('$') > 0):
                in_header = True

            if in_header:
                header.append(line)
            elif after_header:
                data.append(list(map(float, line.split())))

            if line.find('$END') > 0:
                in_header = False
                after_header = True

    # reshape data
    data = np.concatenate([np.array(i) for i in data])
    data = (data[0::2] + 1j * data[1::2]).astype(complex)

    # LCModel-specific conjugation
    if conjugate:
        data = np.conj(data)

    return data, header


#*******************************************#
#   loading LCModel raw fixed header size   #
#*******************************************#
def read_LCModel_raw_hs(filename, header_size=11, conjugate=True):
    """
    Read LCModel raw format with user-specified header size.

    @param filename -- Path to file.
    @param header_size -- Number of header lines.
    @param bool conjugate -- Apply conjugation upon read.

    @returns -- The basis set data/FID and header if possible.
    """
    header = []
    data = []
    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            if i >= header_size: data.append(list(map(float, line.split())))
            else: header.append(line)

    # reshape data
    data = np.concatenate([np.array(i) for i in data])
    data = (data[0::2] + 1j * data[1::2]).astype(complex)

    # LCModel-specific conjugation
    if conjugate:
        data = np.conj(data)

    return data, header


#******************************#
#   load LCModel coord data    #
#******************************#
def read_LCModel_coord(path, coord=True, meta=True):
    """
    Load data based on LCModel coord files.

    @param path -- The path to the files.
    @param coord -- Load concentration estimates.
    @param meta -- Load meta data.

    @returns -- The data.
    """
    metabs, concs, crlbs, tcr = [], [], [], []
    fwhm, snr, shift, phase = None, None, None, None

    # go through file and extract all info
    with open(path, 'r') as file:
        concReader = 0
        miscReader = 0

        for line in file:
            if 'lines in following concentration table' in line:
                concReader = int(line.split(' lines')[0])
            elif concReader > 0:  # read concentration table
                concReader -= 1
                values = line.split()

                # check if in header of table
                if values[0] == 'Conc.':
                    continue
                else:
                    try:  # sometimes the fields are fused together with '+'
                        m = values[3]
                        c = float(values[2])
                    except:
                        if 'E+' in values[2]:  # catch scientific notation
                            c = values[2].split('E+')
                            m = str(c[1].split('+')[1:])
                            c = float(c[0] + 'e+' + c[1].split('+')[0])
                        else:
                            if len(values[2].split('+')) > 1:
                                m = str(values[2].split('+')[1:])
                                c = float(values[2].split('+')[0])
                            elif len(values[2].split('-')) > 1:
                                m = str(values[2].split('-')[1:])
                                c = float(values[2].split('-')[0])
                            else:
                                raise ValueError(f'Could not parse {values}')

                    # append to data
                    metabs.append(m)
                    concs.append(float(values[0]))
                    crlbs.append(int(values[1][:-1]))
                    tcr.append(c)
                    continue

            if 'lines in following misc. output table' in line:
                miscReader = int(line.split(' lines')[0])
            elif miscReader > 0:  # read misc. output table
                miscReader -= 1
                values = line.split()

                # extract info
                if 'FWHM' in values:
                    fwhm = float(values[2])
                    snr = float(values[-1].split('=')[-1])
                elif 'shift' in values:
                    if values[3] == 'ppm':
                        shift = float(values[2][1:])  # negative fuses with '='
                    else:
                        shift = float(values[3])
                elif 'Ph' in values:
                    phase = float(values[1])

    if coord and meta: return metabs, concs, crlbs, tcr, fwhm, snr, shift, phase
    elif coord: return metabs, concs, crlbs, tcr
    elif meta: return fwhm, snr, shift, phase