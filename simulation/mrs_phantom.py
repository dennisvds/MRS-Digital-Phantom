import os
import pandas as pd
import torch
import torch.nn.functional as F
import scipy.ndimage as ndimage

from simulation.signalModel import SignalModel
from simulation.basis import Basis
from utils.loading import load_skeleton_as_subject

class DigitalPhantom:
    """
    Class to generate synthetic MRSI data using a digital phantom skeleton,
    a metabolite basis set, and defined simulation parameters.
    """

    def __init__(self, skeleton, path2metabs, path2basis, simulation_params, gui=None, metabs=[]):
        """
        Initialize the Digital Phantom with skeleton, basis set, and simulation settings.

        Args:
            skeleton (str): Name of the skeleton dataset.
            path2metabs (str): Path to metabolite CSV file.
            path2basis (str): Path to metabolite basis set.
            simulation_params (dict): Dictionary with simulation settings.
            gui (optional): GUI interface object, if used for parameter selection.
            metabs (optional): List of selected metabolites, if not using GUI.
        """
        self.skeleton = skeleton
        self.path2basis = path2basis
        self.gui = gui

        self._log("📥 Initializing Digital Phantom...")

        # Load the skeleton subject and label names
        path2skeleton = os.path.join('./data/skeletons/', self.skeleton)
        self.subject, self.label_names = load_skeleton_as_subject(path2skeleton, self.skeleton)
        self.spacing = self.subject.spacing
        self.affine = self.subject['labels'].affine

        # Store simulation parameters and add voxel spacing info
        self.simulation_params = simulation_params
        self.simulation_params['voxel_spacing'] = self.spacing

        # Load metabolite information
        self.metab_df = pd.read_csv(path2metabs)
        self.metabs = self.gui.basisset_widget.metabolite_checkboxes.get_selected_metabolites() if gui else metabs

        # Create the basis set and signal model
        self._create_basis_and_signal_model()

        # Initialize lipid contamination maps
        self.lipid_label_map = self.init_lipid_labels()
        self.lipid_mask = self.update_lipid_mask(self.lipid_label_map)

    def _log(self, msg):
        """Helper function to log messages to the GUI or stdout."""
        if self.gui:
            self.gui.log_message(msg)
        else:
            print(f"[INFO] {msg}")

    def _create_basis_and_signal_model(self):
        """Helper function to instantiate the basis set and signal model."""
        self.basis_set = Basis(
            self.path2basis, self.metab_df,
            bw=self.simulation_params['bandwidth'],
            points=self.simulation_params['spectral_points'],
            TE=self.simulation_params.get('TE'),
            vendor=self.simulation_params.get('vendor'),
            localization=self.simulation_params.get('localization'),
            metabs=self.metabs
        )

        self.signal_model = SignalModel(
            self.metab_df,
            self.basis_set,
            labels=self.label_names,
            simulation_params=self.simulation_params,
            gui=self.gui
        )

    def get_image_data(self):
        """
        Returns the anatomical image data corresponding to the skeleton.

        Returns:
            np.ndarray: 3D anatomical image data.
        """
        if self.skeleton == 'BigBrainMR':
            return self.subject['T1wImage'].data.numpy()[0]
        elif self.skeleton == 'MRiLab':
            return self.subject['PD'].data.numpy()[0]
        elif self.skeleton == 'Custom':
            return self.subject['T1'].data.numpy()[0]
        else:
            raise NotImplementedError(f"Image type for skeleton '{self.skeleton}' not implemented.")

    def get_label_data(self, one_hot=False):
        """
        Returns the label map for the skeleton.

        Args:
            one_hot (bool): If True, returns a one-hot encoded label map.

        Returns:
            torch.Tensor: Label map (optionally one-hot encoded).
        """
        label_map = self.subject['labels'].data[0]
        if one_hot:
            num_classes = len(torch.unique(label_map))
            return torch.nn.functional.one_hot(label_map.long(), num_classes=num_classes)
        return label_map
    
    def get_B0_data(self):
        """
        Returns the B0 field map data for the skeleton.

        Returns:
            np.ndarray: 3D B0 field map data.
        """
        if 'B0' in self.subject:
            return self.subject['B0'].data.numpy()[0]
        else:
            self._log("⚠️ No B0 field map available for this skeleton.")
            return None

    def simulate_data(self, voi_coords):
        """
        Simulates MRS data for a volume of interest (VOI).

        Args:
            voi_coords (tuple): Coordinates defining the VOI: (x_start, x_end, y_start, y_end, z_start, z_end).

        Returns:
            tuple: Simulated spectrum, component spectra, time points, and ppm scale.
        """
        x0, x1, y0, y1, z0, z1 = voi_coords
        voi_labels = self.get_label_data(one_hot=True)[x0:x1, y0:y1, z0:z1]
        voi_lipid_mask = self.lipid_mask[x0:x1, y0:y1, z0:z1]

        # If B0 map is available, pass it to the signal model
        if 'B0' in self.subject:
            voi_B0_map = self.get_B0_data()[x0:x1, y0:y1, z0:z1]

        spec, components = self.signal_model.simulate_mrs_data(voi_labels, voi_lipid_mask, B0_map=voi_B0_map if 'B0' in self.subject else None)
        return spec, components, self.basis_set.t, self.basis_set.ppm

    def init_lipid_labels(self):
        """
        Initializes the lipid label map based on available brain or fat labels.

        Returns:
            torch.Tensor: Binary mask highlighting lipid regions.
        """
        if 'Fat' in self.label_names:
            self._log("🧬 Generating lipid mask from explicit fat label...")
            fat_idx = self.label_names.index('Fat')
            return (self.get_label_data() == fat_idx).float()

        self._log("⚠️ No explicit fat label found, generating lipid mask from whole-brain mask...")
        whole_brain_mask = self.subject['brain_mask'].data[0]
        eroded = F.max_pool3d(whole_brain_mask.unsqueeze(0).unsqueeze(0), 3, 1, 1)
        eroded = (eroded == 1).float().squeeze(0).squeeze(0)
        return (eroded - whole_brain_mask > 0).float()

    def update_lipid_mask(self, lipid_mask):
        """
        Smooth and normalize the lipid mask to simulate lipid contamination.

        Args:
            lipid_mask (torch.Tensor): Raw lipid label map.

        Returns:
            torch.Tensor: Smoothed and normalized lipid mask.
        """
        sigma = self.simulation_params['lipid_sigma']
        spacing = self.simulation_params['voxel_spacing'][0]
        sigma_voxel = sigma / spacing

        smoothed = ndimage.gaussian_filter(lipid_mask.numpy(), sigma=sigma_voxel, mode='nearest')
        smoothed = torch.tensor(smoothed)
        return smoothed / smoothed.max() if smoothed.max() > 0 else smoothed

    def update_phantom(self, simulation_params):
        """
        Update the phantom with new simulation parameters and update the lipid mask if necessary.

        Args:
            simulation_params (dict): Updated simulation parameters.
        """
        old_value = getattr(self, 'simulation_params', {}).get('lipid_sigma')
        new_value = simulation_params.get('lipid_sigma')
        needs_update = old_value != new_value or not hasattr(self, 'simulation_params')

        if not hasattr(self, 'simulation_params'):
            self._log("🆕 No previous simulation parameters found. Will update lipid mask.")
        elif needs_update:
            self._log(f"🔄 lipid_sigma changed from {old_value} to {new_value}. Updating lipid mask.")
        else:
            self._log(f"➡️ lipid_sigma unchanged ({old_value}). Skipping lipid mask update.")

        # Store new simulation parameters
        self.simulation_params = simulation_params
        self.simulation_params['voxel_spacing'] = self.spacing

        # Recreate basis set and signal model
        self._create_basis_and_signal_model()

        # Update the lipid mask if required
        if self.lipid_mask is not None and needs_update:
            self._log("🔄 Updating lipid mask...")
            self.lipid_mask = self.update_lipid_mask(self.lipid_label_map)
        else:
            self._log("➡️ Lipid mask update skipped.")
