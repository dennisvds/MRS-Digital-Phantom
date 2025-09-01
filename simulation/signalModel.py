####################################################################################################
#                                        simulation.py                                             #
####################################################################################################
# Authors: D.M.J. van de Sande (d.m.j.v.d.sande@tue.nl)
# Created: 03/02/23
# Purpose: Defines methods to generate batches of spectra using the Digital MRS Phantom.
####################################################################################################

#*************#
#   Imports   #
#*************#
import numpy as np
from tqdm import tqdm
import torch
import torch.nn.functional as F
from scipy.integrate import simpson
from scipy.ndimage import gaussian_filter, binary_dilation, distance_transform_edt
import math

# Internal modules
from simulation.macromolecules.simulate_mm import simulate_MMs
from simulation.water.water_sim import generate_sum_fid
from simulation.run_matlab import MatlabRunner
from utils.auxillary import calc_scale_factor, spec2fid, fid2spec
from utils.definitions import H2O_CONCENTRATIONS
from utils.visualization import plot_shim_map_process, plot_spectra

#*************#
# Signal Model #
#*************#
class SignalModel:
    def __init__(self, metab_df, basis, labels=['Background', 'WM', 'GM', 'CSF'], simulation_params=None, gui=None):
        """Initialize the signal model with basis set and simulation parameters."""
        self.metab_df = metab_df
        self.basis = basis
        self.labels = labels
        self.B0 = basis.B0
        self.TE = basis.TE
        self.gui = gui
        self.simulation_params = simulation_params

    def log(self, message: str, color: str = 'black'):
        """Helper function to log messages to GUI or print fallback."""
        if self.gui:
            self.gui.log_message(message, color)
        else:
            print(message)

    def simulate_mrs_data(self, voi_labels, voi_lipid_mask, domain='ppm'):
        """Simulate an MRS spectrum from input VOI labels and masks."""
        if self.simulation_params is None:
            self.log("❌ Simulation parameters not set. Please provide simulation parameters.", 'red')
            raise ValueError("Simulation parameters not set.")

        spec = torch.zeros((self.basis.n), dtype=torch.complex64)
        voi_labels = voi_labels.to(torch.complex64)
        voi_labels_flat = voi_labels.view(-1, voi_labels.shape[-1])

        # Calculate the voxel volume in mm^3
        voxel_volume = torch.prod(torch.tensor(self.simulation_params['voxel_spacing']))

        # ---- Shim map simulation ---- #
        self.log("Simulating shim map...")
        self.shim_map = self.simulate_shim_map(voi_labels)

        # ---- Metabolite signal ---- #
        self.log("Simulating metabolites...")
        metab_spectra = self.simulate_metabs()
        spec_metabs = torch.matmul(voi_labels_flat, metab_spectra)
        spec_metabs = spec_metabs.view(*voi_labels.shape[:3], -1)
        spec_metabs = self.apply_shim_map(spec_metabs, self.shim_map)
        spec_metabs *= voxel_volume  # Scale by voxel volume
        spec_metabs = spec_metabs.sum((0, 1, 2))  # Sum over the VOI dimensions
        spec += spec_metabs

        # ---- Macromolecule signal ---- #
        self.log("Simulating macromolecules...")
        mm_specs = self.simulate_mms()
        spec_mms = torch.matmul(voi_labels_flat, mm_specs)
        spec_mms = spec_mms.view(*voi_labels.shape[:3], -1)
        spec_mms *= self.simulation_params['mm_level'] * 1e3
        spec_mms = self.apply_shim_map(spec_mms, self.shim_map)
        spec_mms *= voxel_volume  # Scale by voxel volume
        spec_mms = spec_mms.sum((0, 1, 2))  # Sum over the VOI dimensions
        spec += spec_mms

        # ---- Lipid signal ---- #
        self.log("Simulating lipids...")
        lipids_spec = self.simulate_lipids(voi_labels, voi_lipid_mask)
        spec_lipids = voi_lipid_mask[..., None] * lipids_spec[None, None, None, :]
        spec_lipids *= self.simulation_params['lipid_amp_factor']
        spec_lipids = self.apply_shim_map(spec_lipids, self.shim_map)
        spec_lipids *= voxel_volume  # Scale by voxel volume
        spec_lipids = spec_lipids.sum((0, 1, 2))  # Sum over the VOI dimensions
        spec += spec_lipids

        # ---- Water signal ---- #
        self.log("Simulating water...")
        specs_water = self.simulate_water()
        spec_water = torch.matmul(voi_labels_flat, specs_water)
        spec_water = spec_water.view(*voi_labels.shape[:3], -1)
        spec_water = self.apply_shim_map(spec_water, self.shim_map)
        spec_water *= voxel_volume  # Scale by voxel volume
        spec_water = spec_water.sum((0, 1, 2))  # Sum over the VOI dimensions
        spec += spec_water

        # ---- Add Noise ---- #
        self.log("Simulating noise...")
        noise_std = self.simulation_params['noise_level'] * 2e4  # Empirical scaling factor to match noise level
        spec_noise = (torch.randn(spec.shape) + 1j * torch.randn(spec.shape)) * noise_std
        spec += spec_noise

        components = {
            'Metabs': spec_metabs,
            'MMs': spec_mms,
            'Lipids': spec_lipids,
            'Water': spec_water,
            'Noise': spec_noise
        }

        self.log("✅ Simulation complete!")
        return spec, components

    def simulate_metabs(self):
        # Remove metabolites from metabolite dataframe that are not in the basis set
        self.log("Removing metabolites from dataframe that are not in the basis set...")
        self.log(f"{len(set(self.metab_df['Metabolite'].unique()) - set(self.basis.names))} metabolites removed: {set(self.metab_df['Metabolite'].unique()) - set(self.basis.names)}", 'red')

        metab_df = self.metab_df[self.metab_df['Metabolite'].isin(self.basis.names)]

        fids = torch.tensor(self.basis.fids[torch.newaxis, :, :])        # (1, N, M)
        t = torch.tensor(self.basis.t[torch.newaxis, :, torch.newaxis])  # (1, N, 1)
        f = torch.tensor(self.basis.f[torch.newaxis, :, torch.newaxis])  # (1, N, 1)

        # # Creatine index for peak area estimation
        # cr_idx = self.basis.names.index('Cr')   
        # # Extract creatine spectrum from basis set
        # cr_fid = fids[0,:,cr_idx]
        # cr_spec = torch.fft.fftshift(torch.fft.fft(cr_fid, dim=0), dim=0)
        # # Find the methyl peak area (-CH3 group with 3 protons, ~ 3.0 ppm)
        # ch3_spec = cr_spec[(self.basis.ppm > 2.9) & (self.basis.ppm < 3.2)]
        # ppm_axis = self.basis.ppm[(self.basis.ppm > 2.9) & (self.basis.ppm < 3.2)]
        # ch3_area = simpson(ch3_spec.real, x=ppm_axis)
        # one_proton_area = ch3_area / 3.0  # Area under the curve for one proton

        # Extract metabolite concentrations for the given labels
        metab_specs = []
        for label in self.labels:
            # Subset dataframe based on tissue label
            df = metab_df[metab_df['Tissue'] == label]

            # If no metabolites are present in the dataframe for this label, set the spectrum to zero
            if df.shape[0] == 0:
                metab_specs.append(torch.zeros((1, self.basis.n), dtype=torch.complex64))
                continue

            # Concentrations
            concs_mean = torch.tensor(df['Conc_mean'].values)
            concs_std = torch.tensor(df['Conc_std'].values)
            concs_mean[torch.isnan(concs_mean)] = 0.0  # Set NaN to zero
            concs_std[torch.isnan(concs_std)] = 0.0    # Set NaN to zero
            concs = torch.normal(concs_mean, concs_std)  # Sample from normal distribution

            # T2 relaxation (Lorentzian linewidth)
            t2s = torch.tensor(df['T2'].values)
            gamma = 1/(torch.pi * t2s * 1e-3)
            # Gaussian linewidth
            sigma = torch.ones_like(concs) * 0.0
            # Frequency shift per metabolite
            eps = torch.ones_like(concs) * 0.0

            # set nan values to zero
            gamma[torch.isnan(gamma)] = 0

            # Reshape
            fid = fids * torch.exp(-1j * eps * t - gamma * t - (sigma*t)**2)
            spec = torch.fft.fftshift(torch.fft.fft(fid, dim=1), dim=1)
            
            # Multiply with concentrations
            specs = concs * spec
            specs = specs.sum(dim=-1)

            # Add to list
            metab_specs.append(specs)
        

        # Convert to tensor (torch.complex64)
        metab_specs = torch.stack(metab_specs, dim=0).squeeze()
        metab_specs = metab_specs.to(torch.complex64)

        return metab_specs


    def simulate_mms(self):
        mm_spec = simulate_MMs(tesla=self.B0, 
                                           basis=self.basis,
                                           TE=self.TE, 
                                           TR=self.simulation_params["TR"], 
                                           json_file=self.simulation_params['mm_json_file'],
                                           )
        
        mm_specs = torch.zeros((len(self.labels), self.basis.n), dtype=torch.complex64)
        mm_specs[1, :] = mm_spec[0, :] # WM
        mm_specs[2, :] = mm_spec[1, :] # GM

        return mm_specs

    def simulate_lipids(self, voi_labels, voi_lipid_mask):

        matlab_runner = MatlabRunner()
        lipid_fid = matlab_runner.generate_lipid_basis(basis=self.basis, linewidth=5.0)
        lipid_components = lipid_fid.shape[0]

        # Time vector
        t = torch.tensor(self.basis.t, dtype=torch.float64)

        # Set linewidth and phase parameters
        linewidth_min, linewidth_max = self.simulation_params['lipid_lw_min'], self.simulation_params['lipid_lw_max']
        phase_min, phase_max = self.simulation_params['lipid_phase_min'], self.simulation_params['lipid_phase_max']

        # Sample random values from uniform distributions
        linewidths = torch.empty(lipid_components).uniform_(linewidth_min, linewidth_max)  # (components,)
        phases = torch.empty(lipid_components).uniform_(phase_min, phase_max)    # (components,)    

        # Expand dimensions for broadcasting
        linewidths = linewidths[:, None]  # (components, 1)
        phases = phases[:, None]          # (components, 1)

        # Compute per-component decay and phase
        decay = torch.exp(-np.pi * linewidths * t)  # (components, time)
        phase_shift = torch.exp(1j * phases * np.pi / 180.0)  # (components, 1), radians

        # Apply decay and phase separately to each component
        lipid_fid_proc = lipid_fid * decay * phase_shift

        # Convert to spectrum
        lipid_spec = fid2spec(lipid_fid_proc, axis=-1)

        # Sum over the components
        lipid_spec = lipid_spec.sum(dim=0)  # (N,)


        return lipid_spec
    
    def simulate_water(self):
        # Simulation parameters
        WF = self.simulation_params['water_amp_factor']*10e-3  # Water scaling factor
        n_comp = 5
        amplitudes = WF * np.array([30, 80, 180, 70, 40])  # amplitudes
        freq_ppm    = np.array([4.90, 4.82, 4.70, 4.62, 4.50])
        
        damping_hz  = np.random.uniform(self.simulation_params['water_damping_min'], 
                                        self.simulation_params['water_damping_max'], n_comp)
        
        phase_deg   = np.random.uniform(self.simulation_params['water_phase_min'], 
                                       self.simulation_params['water_phase_max'], n_comp)

        water_fid_lorentz = generate_sum_fid(self.basis.t, amplitudes, 
                                             damping_hz, 
                                             freq_ppm, 
                                             phase_deg, 
                                             ref_freq_hz=self.basis.cf, 
                                             decay_type='lorentzian', scaling=WF)
        
        water_fid_gauss = generate_sum_fid(self.basis.t, amplitudes, 
                                             damping_hz, 
                                             freq_ppm, 
                                             phase_deg, 
                                             ref_freq_hz=self.basis.cf, 
                                             decay_type='gaussian', scaling=WF)
        
        # Combine the two FIDs
        water_fid = water_fid_lorentz / 2 + water_fid_gauss / 2
        # Convert to complex64
        water_fid = torch.tensor(water_fid, dtype=torch.complex64)
        # Transform to frequency domain
        water_spec = fid2spec(water_fid, axis=-1)  # Transform to frequency domain
        # Create water spec for each label
        water_specs = []
        for label in self.labels:
            conc = H2O_CONCENTRATIONS[label]
            # Scale the water spectrum by the concentration of the labels
            water_spec_label = water_spec * conc
            water_specs.append(water_spec_label)
        # Stack the water spectra for each label
        water_specs = torch.stack(water_specs, dim=0)  # Shape: (l, x, y, z, n_points)
        water_specs = water_specs.to(torch.complex64)  # Convert to complex64

        return water_specs
    
    def simulate_shim_map(self, voi_labels, seed=None):
        """
        Simulate a 3D shim imperfection map with a smooth boundary effect.
        The boundary effect gradually transitions inside and outside the brain.
        """
        # --- unpack simulation parameters ---
        sigma_hz        = self.simulation_params['shim_amplitude_hz']
        corr_length_mm  = self.simulation_params['shim_corr_length']
        boundary_amp    = self.simulation_params['shim_boundary_amp_factor']
        boundary_smooth = self.simulation_params['shim_boundary_smoothing']
        voxel_spacing   = self.simulation_params['voxel_spacing'][0]

        if seed is not None:
            np.random.seed(seed)

        # Ensure real-valued labels
        voi_labels = np.array(voi_labels.real)

        # Step 1: collapse one-hot to int label map
        label_map = np.argmax(voi_labels, axis=3)  # shape (nx,ny,nz)
        nx, ny, nz = label_map.shape

        # Step 2: base smooth random field
        noise   = np.random.randn(nx, ny, nz)
        sigma_v = corr_length_mm / voxel_spacing
        base    = gaussian_filter(noise, sigma=sigma_v, mode='reflect')
        base   /= base.std()
        base   *= sigma_hz

        # Step 3: masks
        background_mask = (label_map == 0)
        tissue_mask     = (label_map > 0)

        # Step 4: signed distance transform
        distance_to_boundary = (distance_transform_edt(background_mask, sampling=voxel_spacing) - 
                                distance_transform_edt(tissue_mask,sampling=voxel_spacing))

        # Step 5: smooth boundary transition — no hard cutoff
        boundary_weight = np.exp(-np.abs(distance_to_boundary) / boundary_smooth)
        boundary_ampl   = boundary_amp * boundary_weight

        # Step 6: apply modulation
        final_map = base + boundary_ampl * base

        # Step 7: zero outside brain
        final_map[background_mask] = 0.0

        # Optional plot
        # plot_shim_map_process(label_map, base, boundary_ampl, final_map, show=True)

        return torch.tensor(final_map)

    
    def apply_shim_map(self, specs, shim_map):
        """
        Apply the shim imperfection map to the spectra.
        
        Args:
            specs (torch.Tensor): The input spectra in shape (x,y,z,n_points).
            shim_map (torch.Tensor): The shim imperfection map in shape (x,y,z).
        
        Returns:
            torch.Tensor: The spectra after applying the shim imperfection map.
        """
        # Reshape shim_map to match the shape of specs
        shim_map = shim_map.unsqueeze(-1)

        # Transform specs to time domain
        fids = spec2fid(specs, axis=-1)
        t_axis = self.basis.t

        phase_modulation = torch.exp(-1j * 2 * np.pi * shim_map * t_axis)
        fids_mod = fids * phase_modulation

        # Transform back to frequency domain
        specs_mod = fid2spec(fids_mod, axis=-1)
        return specs_mod






