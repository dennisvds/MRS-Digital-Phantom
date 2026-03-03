####################################################################################################
#                                             basis.py                                             #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#          D.M.J. van de Sande (d.m.j.v.d.sande@tue.nl)                                            #
#                                                                                                  #
# Created: 03/02/23                                                                                #
#                                                                                                  #
# Purpose: Defines a main structure for MRS metabolite basis sets. Encapsulates all information    #
#          and holds definitions to compute various aspects.                                       #
#                                                                                                  #
####################################################################################################


#*************#
#   imports   #
#*************#
import numpy as np
import os
import warnings
import re

# own
from fsl_mrs.utils.mrs_io.main import read_basis
from utils.definitions import METABS, METABS_TRANSLATION, GYRO_MAG_RATIO


#**************************************************************************************************#
#                                           Class Basis                                            #
#**************************************************************************************************#
#                                                                                                  #
# The main structure for MRS metabolite basis sets.                                                #
#                                                                                                  #
#**************************************************************************************************#
class Basis():

    #*************************#
    #   initialize instance   #
    #*************************#
    def __init__(self, path2basis, metab_df=None, bw=None, points=None, metabs=[], vendor=None, localization=None, TE=None):
        """
        Main initialization for the Basis class.
        """
        self.path2basis = path2basis
        self.metab_df = metab_df
        self.bw = bw
        self.points = points
        self.metabs = metabs
        self.vendor = vendor
        self.localization = localization
        self.TE = TE
        self.load_and_format_basis()

    def load_and_format_basis(self):
        """
        Loads and formats the basis set according to the current parameters.
        """
        # Load basis set
        basis = read_basis(self.path2basis)
        
        # Rename metabolites
        basis = self.rename_metabolites(basis)
        
        # Define metabolites to ignore
        ignore = self.define_ignore(basis, self.metabs, self.metab_df)
        
        # Set bandwidth and points if not provided
        if self.bw is None and self.points is None:
            self.bw = basis.original_bw
            self.points = basis.original_points
        

        # Reformat the basis set
        basis = self.reformat(basis, self.bw, self.points, ignore=ignore)

        # Extract water FID and ref0ppm (if present)
        if 'H2O' in basis.names:
            water_idx = basis.names.index('H2O')
            # self.water_fid = basis.original_basis_array[:, water_idx]
            # Remove water from the basis set
            basis.remove_fid_from_basis('H2O')
        if 'Ref0ppm' in basis.names:
            ref0ppm_idx = basis.names.index('Ref0ppm')
            # self.ref0ppm_fid = basis.original_basis_array[:, ref0ppm_idx]
            # Remove ref0ppm from the basis set
            basis.remove_fid_from_basis('Ref0ppm')
        
        # Set attributes
        self.basisFSL = basis
        self.names = [n.split('.')[0] for n in basis.names]
        self.n_metabs = len(self.names)
        self.fids = basis._raw_fids
        self.dwelltime = float(basis.original_dwell)
        self.n = basis.original_points
        self.t = np.arange(self.dwelltime, self.dwelltime * (self.n + 1), self.dwelltime)
        self.f = np.arange(-self.bw / 2, self.bw / 2, self.bw / self.n)
        self.ppm = basis.original_ppm_shift_axis
        self.cf = float(basis.cf)
        self.B0 = self.cf / GYRO_MAG_RATIO['1H'] # Currently hardcoded for 1H

        # # Extract vendor, localization, and TE from the filename
        # filename = os.path.basename(self.path2basis)
        # match = re.match(r'LCModel_(.+?)_UnEdited_(.+?)_TE(\d+(?:\.\d+)?)\.BASIS', filename)
        # if match:
        #     self.vendor = match.group(1)
        #     self.localization = match.group(2)
        #     self.TE = float(match.group(3))
        # else:
        #     raise ValueError(f"Filename does not match expected pattern: {filename}")

    def reset_basis(self, metabs=None, bw=None, points=None):
        """
        Reset the basis set with new parameters.
        """
        if metabs is not None:
            self.metabs = metabs
        if bw is not None:
            self.bw = bw
        if points is not None:
            self.points = points
        
        # Reload and reformat the basis set
        self.load_and_format_basis()
    
    #*************************#
    #   get formatted basis   #
    #*************************#
    def reformat(self, basis, bw, points, ignore=[]):
        """
        Reformat the basis set to a given bandwidth and number of points.
        """
        basis._raw_fids = basis.get_formatted_basis(bw, points, ignore=ignore)
        basis._raw_fids /= np.mean(np.abs(basis._raw_fids))
        basis._dt = 1. / bw
        basis._names = basis.get_formatted_names(ignore=ignore)
        return basis
    
    def rename_metabolites(self, basis):
        for name in basis.names[:]:
            if name not in METABS:
                if name in METABS_TRANSLATION:
                    basis.names[basis.names.index(name)] = METABS_TRANSLATION[name]
                else:
                    print(f"Warning: '{name}' not in the list of translation. Removing from basisset...")
                    basis.remove_fid_from_basis(name)
        
        sorted_indices = sorted(range(len(basis.names)), key=lambda k: basis.names[k].lower())
        new_names = [basis.names[i] for i in sorted_indices]
        new_arrays = basis.original_basis_array[:, sorted_indices]
        for i in range(len(new_names)):
            basis.names[i] = new_names[i]
            basis.original_basis_array[:, i] = new_arrays[:, i]
        
        return basis
    
    def define_ignore(self, basis, metabs, metab_df):
        all_basis_metabs = basis.names

        for m in metabs:
            if m not in all_basis_metabs:
                warnings.warn(f"Metabolite '{m}' not in basisset. Ignoring...")

        ignore = [] if len(metabs) == 0 else [m for m in all_basis_metabs if m not in metabs]

        if metab_df is not None:
            all_df_metabs = metab_df['Metabolite'].unique()
            ignore += [m for m in all_basis_metabs if m not in all_df_metabs]

        # Remove duplicates
        ignore = list(set(ignore))

        # Make sure 'H2O' is not ignored
        if 'H2O' in ignore:
            ignore.remove('H2O')

        print(f"Metabolites to ignore from basisset: {ignore}")
        return ignore
    
    def summary(self):
        print(f"Summary of the basis set:")
        print(f"{self.basisFSL}")
        print(f"Loaded from: {self.path2basis}")
