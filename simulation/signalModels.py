####################################################################################################
#                                             signalModels.py                                      #
####################################################################################################
#                                                                                                  #
# Authors: J. P. Merkofer (j.p.merkofer@tue.nl)                                                    #
#          D.M.J. van de Sande (d.m.j.v.d.sande@tue.nl)                                            #
#                                                                                                  #
# Created: 03/02/23                                                                                #
#                                                                                                  #
# Purpose: Defines MRS signal models to generate spectra.                                          #
#                                                                                                  #
####################################################################################################

import numpy as np
from scipy.fft import ifft, ifftshift, fftshift, fft

# Own
from simulation.simulate_mms import simulate_MMs
from utils.definitions import lipid_params
from utils.definitions import residual_water_params as res_params
from utils.auxillary import randomWalk, randomPeak, baseline_init

class SignalModel():
    def __init__(self, basis, gauss_broadening=10, freq_shift=0, mm_scale=1000, noise_std=0.1, res_water_scale=10, lipid_scale=10, baseline_order=2, baseline_scale=400):
        # Extract parameters from basis
        self.fids = basis.fids
        self.t = basis.t
        self.f = basis.f
        self.n = basis.n
        self.cf = basis.cf
        self.ppm = basis.ppm
        self.basis = basis

        # User inputs
        self.gauss_broadening = gauss_broadening
        self.freq_shift = freq_shift
        self.mm_scale = mm_scale
        self.noise_std = noise_std
        self.res_water_scale = res_water_scale
        self.lipid_scale = lipid_scale
        self.baseline_order = baseline_order
        self.baseline_scale = baseline_scale

        # Compute water residual
        self.res_water_spec = randomWalk(self.ppm, rw_range=res_params['ppm_range'], scale=res_params['scale'], 
                                         smooth=res_params['smooth'], ylim=res_params['ylim'])
        
        # Compute baseline
        if self.baseline_order is not None:
            self.baseline_spec = baseline_init(points=self.n, order=self.baseline_order, first=0, last=self.n)
            self.baseline_spec = np.sum(self.baseline_spec, axis=-1) * self.baseline_scale


    def forward(self, phantom_data, labels, lipid_mask=None):
        # Prepare parameters
        params = self.prepare_params(phantom_data)

        # Split parameters
        con = params[0][:, np.newaxis, :]       # Shape: (bs, num_metabs, 1)
        eps = params[1][:, np.newaxis, :]       # Shape: (bs, num_metabs, 1)
        gamma = params[2][:, np.newaxis, :]     # Shape: (bs, num_metabs, 1)
        sigma = params[3][:, np.newaxis, :]     # Shape: (bs, num_metabs, 1)

        # Set NaN values to zero to avoid errors 
        # This will only simulate a metabolite when all parameters are available 
        nan_mask = np.isnan(con) | np.isnan(eps) | np.isnan(gamma) | np.isnan(sigma)
        con[nan_mask] = 0
        eps[nan_mask] = 0
        gamma[nan_mask] = 0
        sigma[nan_mask] = 0

        # Reshape FIDs, t, and f
        fids = self.fids[np.newaxis, :, :]      # Shape: (1, n_points, num_metabs)
        t = self.t[np.newaxis, :, np.newaxis]   # Shape: (1, n, 1)
        f = self.f[np.newaxis, :, np.newaxis]   # Shape: (1, 2048, 1)

        # Compute: M_m(v; a_m, w_m) = FFT(m_m(t)* exp(-i*w_m*t - a_m*t - g*t^2))
        fids = self.fids * np.exp(-1j * eps * t - gamma * t - (sigma * t) ** 2)
        specs = fftshift(fft(fids, axis=1), axes=1)

        # Multiply metabolite spectra by concentrations and sum over metabolites
        specs = con * specs
        specs = np.sum(specs, axis=-1)

        # Add macromolecular baseline
        mm_specs = simulate_MMs(tesla=3.0, basis=self.basis, labels=labels, scaling=self.mm_scale)
        specs = specs + mm_specs

        # Add water residual
        res_water_spec = np.broadcast_to(self.res_water_spec, (phantom_data.shape[0], self.n))
        specs = specs + res_water_spec * self.res_water_scale

        # Add lipid contamination
        if lipid_mask is not None:
            lipid_idxs = [np.argmin(np.abs(self.ppm - lipid_params['ppm_pos'][i])) for i in range(len(lipid_params['ppm_pos']))]
            lipid_widths = lipid_params['ppm_widths'] * self.cf
            peaks = [randomPeak(waveLength=self.n, amp=lipid_params['amps'][i], pos=lipid_idxs[i], width=lipid_widths[i], phase=lipid_params['phases'][i]) 
                    for i in range(len(lipid_params['amps'])-1)]
            lipid_specs = np.sum(peaks, axis=0)
            lipid_specs = np.broadcast_to(lipid_specs, (phantom_data.shape[0], self.n))
            lipid_specs = lipid_specs * lipid_mask[:, np.newaxis]
            specs = specs + lipid_specs * self.lipid_scale

        # Add baseline
        if self.baseline_order is not None:
            specs = specs + self.baseline_spec

        # Add noise
        noise_real = np.random.normal(0, self.noise_std, specs.shape)
        noise_imag = np.random.normal(0, self.noise_std, specs.shape)
        noise = noise_real + 1j * noise_imag
        specs = specs + noise

        return specs
    
    def prepare_params(self, phantom_data):
        # Concentrations
        con = phantom_data[..., 1]

        # Relaxation times
        t1_mean = phantom_data[..., 2]      # Currently not used in the simulation framework
        t2_mean = phantom_data[..., 3]

        # Lorentzian decay constant based on T2
        gamma = 1/(np.pi * t2_mean * 1e-3)
        # Gaussian decay constant based on user input
        sigma = np.ones_like(gamma) * self.gauss_broadening

        # Frequency shift per metabolite
        eps = np.ones_like(gamma) * self.freq_shift

        # Pack parameters
        params = [con, eps, gamma, sigma]

        return params