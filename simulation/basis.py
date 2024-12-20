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
import warnings

# own
from loading.loadBasis import loadBasisAsFSL
from utils.definitions import METABS, METABS_TRANSLATION


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
    def __init__(self, path2basis, fmt=None, bw=2000, points=2048, metabs=[], phantom_metabs=[]):
        """
        Main initialization for the Basis class.

        Args:
        - path2basis (str): The path to the basis set folder.
        - fmt (str): The format of the basis set.
        - bw (int): The bandwidth to reformat to.
        - points (int): The number of points to reformat to.
        - metabs (list): The metabolites to ignore from the basis set. If empty, all metabolites are used.
        - phantom_metabs (list): The metabolites to ignore from the basis set based on the phantom. If empty, all metabolites are used.
        """
        # Load basis set
        basis = loadBasisAsFSL(path2basis, fmt=fmt)

        # Rename metabolites
        basis = self.rename_metabolites(basis)

        # Define metabolites to ignore
        ignore = self.define_ignore(basis, metabs, phantom_metabs)


        # Reformat basis set
        basis = self.reformat(basis, bw, points, ignore=ignore)

        # Set attributes
        self.basisFSL = basis
        self.names = [n.split('.')[0] for n in basis.names]
        self.n_metabs = len(self.names)
        self.fids = basis._raw_fids
        self.bw = basis.original_bw
        self.dwelltime = float(basis.original_dwell)
        self.n = basis.original_points
        self.t = np.arange(self.dwelltime, self.dwelltime * (self.n + 1), self.dwelltime)
        self.f = np.arange(- self.bw / 2, self.bw / 2, self.bw / self.n)
        self.ppm = basis.original_ppm_shift_axis
        self.cf = float(basis.cf)


    #*************************#
    #   get formatted basis   #
    #*************************#
    def reformat(self, basis, bw, points, ignore=[]):
        """
        Reformat the basis set to a given bandwidth and number of points.

        @param basis -- The basis set to reformat.
        @param bw -- The bandwidth to reformat to.
        @param points -- The number of points to reformat to.
        @param ignore -- The metabolites to ignore.

        @returns -- The reformatted basis set.
        """
        basis._raw_fids = basis.get_formatted_basis(bw, points, ignore=ignore)
        basis._raw_fids /= np.mean(np.abs(basis._raw_fids))
        basis._dt = 1. / bw
        basis._names = basis.get_formatted_names(ignore=ignore)
        return basis
    
    def rename_metabolites(self, basis):
        # Translation of metabolite names
        for name in basis.names:
            if name not in METABS:
                # Translate the name
                if name in METABS_TRANSLATION:
                    basis.names[basis.names.index(name)] = METABS_TRANSLATION[name]
                else:
                    print(f"Warning: '{name}' not in the list of translation. Removing from basisset...")
                    basis.remove_fid_from_basis(name)

        # Sort metabolites alphabetically
        sorted_indices = sorted(range(len(basis.names)), key=lambda k: basis.names[k].lower())
        new_names = [basis.names[i] for i in sorted_indices]
        new_arrays = basis.original_basis_array[:, sorted_indices]
        for i in range(len(new_names)):
            basis.names[i] = new_names[i]
            basis.original_basis_array[:, i] = new_arrays[:, i]

        return basis
    
    def define_ignore(self, basis, metabs, phantom_metabs):
        all_basis_metabs = basis.names

        # Indicate when a metabolite is not in the basisset
        for m in metabs:
            if m not in all_basis_metabs:
                warnings.warn(f"Metabolite '{m}' not in basisset. Ignoring...")
            
        # Select metabolites to ignore
        if len(metabs) == 0:
            ignore = []
        else:
            ignore = [m for m in all_basis_metabs if m not in metabs]
            print(f"Metabolites to ignore from basisset: {ignore}")

        # Add metabolites to ignore depening on what is present in the phantom
        if len(phantom_metabs) > 0:
            ignore += [m for m in all_basis_metabs if m not in phantom_metabs]


        return ignore

