import json
import numpy as np
import torch
from scipy.signal import hilbert
from scipy.special import voigt_profile
from scipy.stats import norm, cauchy
from simulation.macromolecules.bloch_simulation import mm_bloch_simulation

#%% Helper function to load parameters from a JSON file
def load_mm_params(json_file, tesla):
    with open(json_file, "r") as f:
        params = json.load(f)
    
    # Calculate B0 in Hz from tesla
    B0 = params["B0_factor"] * tesla
    
    # Convert positions (in ppm-like factors) to Hz and apply offset
    params["position"] = [p * B0 - (params["ppm_offset"] * B0) for p in params["ppm_position"]]
    
    # Convert T1 and T2 values from ms to seconds
    params["T1_gm"] = [t / 1000.0 for t in params["T1_gm"]]
    params["T1_wm"] = [t / 1000.0 for t in params["T1_wm"]]
    params["T2_gm"] = [t / 1000.0 for t in params["T2_gm"]]
    params["T2_wm"] = [t / 1000.0 for t in params["T2_wm"]]
    
    return params

#%% Macromolecule class for defining macromolecules
class Macromolecule:
    """
    Macromolecule characteristics to be used in attenuation calculations.
    """
    def __init__(self, T1_pure_GM, T1_pure_WM, T2_pure_GM, T2_pure_WM,
                 lorentz_width_gm, lorentz_width_wm, gauss_width_gm, gauss_width_wm,
                 Bloch_magnetization_gm, Bloch_magnetization_wm, pos, name, scale, J_coupling_1, J_coupling_2):
        self.T1_GM = T1_pure_GM
        self.T1_WM = T1_pure_WM
        self.T2_GM = T2_pure_GM
        self.T2_WM = T2_pure_WM
        self.pos = pos
        self.name = name
        self.scale = scale  # Relative to MM09
        self.lorentz_width_gm = lorentz_width_gm
        self.lorentz_width_wm = lorentz_width_wm
        self.gauss_width_gm = gauss_width_gm
        self.gauss_width_wm = gauss_width_wm
        self.Bloch_magnetization_gm = Bloch_magnetization_gm
        self.Bloch_magnetization_wm = Bloch_magnetization_wm
        self.J1 = J_coupling_1
        self.J2 = J_coupling_2
        self.T1 = None  # to be set per tissue
        self.T2 = None  # to be set per tissue

    def set_T1(self, tissue):
        if tissue == 'WM':
            self.T1 = self.T1_WM
        elif tissue == 'GM':
            self.T1 = self.T1_GM
        else:
            raise ValueError('Tissue not recognized')

    def set_T2_and_linewidths(self, tissue):
        if tissue == 'WM':
            self.T2 = self.T2_WM
            self.lorentz_width = self.lorentz_width_wm
            self.gauss_width = self.gauss_width_wm
        elif tissue == 'GM':
            self.T2 = self.T2_GM
            self.lorentz_width = self.lorentz_width_gm
            self.gauss_width = self.gauss_width_gm
        else:
            raise ValueError('Tissue not recognized')
        
    def set_bloch_magnetization(self, tissue):
        if tissue == 'WM':
            self.Bloch_magnetization = self.Bloch_magnetization_wm
        elif tissue == 'GM':
            self.Bloch_magnetization = self.Bloch_magnetization_gm
        else:
            raise ValueError('Tissue not recognized')

#%% Voigt lines function (same as before)
def gaussian(pos, fwhm, x):
    """
    Creates Gaussian lines for the MM model.
    Parameters
    ----------
    pos : float
        Position of the peak in number of points.
    fwhm : float
        Full width at half maximum (FWHM) of the Gaussian.  
    x : array
        number of points.
    Returns
    -------
    gauss : array
        Gaussian profile line.
    """
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))  # Convert FWHM to standard deviation
    gauss = norm.pdf(x, loc=pos, scale=sigma)
    return gauss

def lorentzian(pos, fwhm, x):
    """
    Creates Lorentzian lines for the MM model.
    Parameters
    ----------
    pos : float
        Position of the peak in number of points.
    fwhm : float
        Full width at half maximum (FWHM) of the Lorentzian.
    x : array
        number of points.
    Returns
    -------
    lorentz : array
        Lorentzian profile line.
    """
   
    gamma = fwhm / 2.0  # Convert FWHM to half-width at half-maximum (HWHM)
    lorentz = cauchy.pdf(x, loc=pos, scale=gamma)
    return lorentz

def voigt(pos, fwhm_lorentz, fwhm_gaussian, x):
    """
    Creates Voigt lines for the MM model.
    Parameters
    ----------
    pos : float
        Position of the peak in number of points.
    fwhm_lorentz : float
        Lorentzian width (FWHM).
    fwhm_gaussian : float
        Gaussian width (FWHM).
    x : array
        number of points.

    Returns
    -------
    voigt_line : array
        Voigt profile line.

    """
    sigma = fwhm_gaussian / (2 * np.sqrt(2 * np.log(2)))
    gamma = fwhm_lorentz / 2.0
    # Use scipy's voigt_profile for accurate Voigt profile calculation
    voigt_line = voigt_profile(x - pos, sigma, gamma)

    return voigt_line



#%% defining sequences to be used for attenuation
class Sequence():
    """
    Sequences where each sequence has its own excitation scheme. The sequences
    are used for attenuating signal in order to simulate the macromolecular 
    baseline.
    """
    
    def __init__(self, name):
        self.name = name
        
    def DIR_semiLASER(TI1, TI2, TE, TR, peak):#  T1, T2, J1, J2):
        T1 = peak.T1
        T2 = peak.T2
        J1 = peak.J1
        J2 = peak.J2
        TE_AHP = TE#(TE/2 + 0.5*(0.000211) + 0.0012*0.25)/2.0 #de Graaf 6.2.7 section of book talks about 
        
        if J2==0:
            j_evolution =  np.cos(np.pi * J1 * TE)
        else:
            j_evolution =  0.5*np.cos(np.pi * J1 * TE) + 0.5*np.cos(np.pi*J2*TE)
#        0.5*np.cos(np.pi*17.5*0.040)
#        0.5*np.cos(np.pi*8*0.040)
#        
#        TE = 40/1000.0
#        beta_3 = 0.5*(0.5*np.cos(np.pi*(-17.5)*TE) + 0.5*np.cos(np.pi*9*TE)) #J-beta-3 2.663ppm
#        beta_2 = 0.5*(0.5*np.cos(np.pi*(-17.5)*TE) + 0.5*np.cos(np.pi*5.5*TE)) #J-beta-2 2.695ppm
#        print(beta_3, beta_2)
        signal = (1-2*np.exp(-TI2/T1)+2*np.exp(-TI1/T1)*np.exp(-TI2/T1)+np.exp(-TR/T1))*np.exp(-TE_AHP/T2)*j_evolution
        return signal
    ###################################
    def IR_semiLASER(TI, TE, TR, peak):
        T1 = peak.T1
        T2 = peak.T2
        J1 = peak.J1
        J2 = peak.J2

        TE = TE#/4.0
        if J2==0:
            j_evolution =  np.cos(np.pi * J1 * TE)
        else:
            j_evolution =  0.5*np.cos(np.pi * J1 * TE) + 0.5*np.cos(np.pi*J2*TE)

        signal = (1-2*np.exp(-TI/T1)+np.exp(-TR/T1))*np.exp(-TE/T2)*j_evolution
        return signal
    ###################################
    def semiLASER(TE, TR, peak):
        T1 = peak.T1
        T2 = peak.T2

        J1 = peak.J1
        J2 = peak.J2
        
        TE = TE #TE/4.0 
        if J2==0:
            j_evolution =  np.cos(np.pi * J1 * TE)
        else:
            j_evolution =  0.5*np.cos(np.pi * J1 * TE) + 0.5*np.cos(np.pi*J2*TE)
        signal = (1+np.exp(-TR/T1))*np.exp(-TE/T2)*j_evolution
        return signal
    ####################################
    def PRESS(TE, TR, peak):
        T1 = peak.T1
        T2 = peak.T2

        J1 = peak.J1
        J2 = peak.J2
        
        TE = TE #TE/4.0 
        if J2==0:
            j_evolution =  np.cos(np.pi * J1 * TE)
        else:
            j_evolution =  0.5*np.cos(np.pi * J1 * TE) + 0.5*np.cos(np.pi*J2*TE)
        signal = (1+np.exp(-TR/T1))*np.exp(-TE/T2)*j_evolution
        return signal
    ###################################
    def IR_STEAM(TI, TE, TM, TR, peak):
        T1 = peak.T1
        T2 = peak.T2
        J1 = peak.J1
        J2 = peak.J2
        if J2==0:
            j_evolution =  np.cos(np.pi * J1 * TE)
        else:
            j_evolution =  0.5*np.cos(np.pi * J1 * TE) + 0.5*np.cos(np.pi*J2*TE)
        
        signal = ((1-2*np.exp(-TI/T1)+np.exp(-TR/T1))* #inversion portion
                 (0.5*1*np.exp(-TM/T1))*np.exp(-TE/T2)*j_evolution) # the 1 is for the sin^3 dependence in the signal equation
        return signal
    ###################################
    def STEAM(TE, TM, TR, peak):
        T1 = peak.T1
        T2 = peak.T2#/4.0
        J1 = peak.J1
        J2 = peak.J2
        if J2==0:
            j_evolution =  np.cos(np.pi * J1 * TE)
        else:
            j_evolution =  0.5*np.cos(np.pi * J1 * TE) + 0.5*np.cos(np.pi*J2*TE)
        signal = (0.5*1*np.exp(-TM/T1))*np.exp(-TE/T2)*(1-np.exp(-TR/T1))*j_evolution
        return signal
    ###################################
    def FID(TE, TR, peak):
        T2_star = peak.T2
        T1 = peak.T1
        
        signal = np.exp(-TE/T2_star)*(1 + np.exp(-TR/T1))
        return signal

# %% Voigt linewidth approximation
# Approximation based on: https://www.sciencedirect.com/science/article/pii/0022407377901613?via%3Dihub (0.02% accuracy)
# f_V ≈ 0.5343 f_L + sqrt(0.2169 f_L^2 + f_G^2)

def calc_gaussian_width(f_voigt, f_lorentz, clip=True):
    """
    Calculate the Gaussian FWHM (f_G) from Voigt and Lorentzian FWHMs.
    """
    term = (f_voigt - 0.5343 * f_lorentz)**2 - 0.2169 * f_lorentz**2
    if np.any(term < 0):
        if clip:
            print("Warning: Negative term in Gaussian width calculation. Clipping to avoid complex linewidths.")
            term = np.clip(term, a_min=1, a_max=None)
        else:
            raise ValueError("Invalid combination: sqrt of negative number. Use clip=True to avoid this.")
    return np.sqrt(term)

def calc_lorentzian_width(f_voigt, f_gauss, clip=True):
    """
    Calculate the Lorentzian FWHM (f_L) from Voigt and Gaussian FWHMs.
    """
    A = 0.0686
    B = -1.0686 * f_voigt
    C = f_voigt**2 - f_gauss**2
    discriminant = B**2 - 4 * A * C
    if np.any(discriminant < 0):
        if clip:
            print("Warning: Negative discriminant in Lorentzian width calculation. Clipping to avoid complex linewidths.")
            discriminant = np.clip(discriminant, a_min=1, a_max=None)
        else:
            raise ValueError("Negative discriminant: no real solution for f_L.")
    return (-B + np.sqrt(discriminant)) / (2 * A)

def calc_voigt_width(f_lorentz, f_gauss):
    """
    Calculate the Voigt FWHM (f_V) from Lorentzian and Gaussian FWHMs.
    """
    return 0.5343 * f_lorentz + np.sqrt(0.2169 * f_lorentz**2 + f_gauss**2)

#%% Compute linewidth components
def calc_lw_components(mm_linewidths, MM_T2s, water_linewidth, water_T2):
    mm_lorentz_widths = 1 / (np.pi * np.array(MM_T2s))
    water_lorentz_width = 1 / (np.pi * np.array(water_T2))

    water_gauss_width = calc_gaussian_width(water_linewidth, water_lorentz_width)

    mm_gauss_widths = calc_gaussian_width(mm_linewidths, mm_lorentz_widths)

    mm_overlap_sq = mm_gauss_widths**2 - water_gauss_width**2
    if np.any(mm_overlap_sq < 0):
        print("Warning: Negative values detected in quadratic subtraction for Gaussian overlap widths. Clipping to avoid complex linewidths.")
        mm_overlap_sq = np.clip(mm_overlap_sq, a_min=1, a_max=None)

    mm_overlap_widths = np.sqrt(mm_overlap_sq)

    lorentz_width = mm_lorentz_widths
    gauss_component = mm_overlap_widths

    # overlap_component = np.array(linewidths) - lorentz_width - np.array(shim_components)
    # gauss_component = np.array(shim_components) + overlap_component
    # gauss_component = np.clip(gauss_component, a_min=1, a_max=None) # Ensure positive values with at least 1 Hz
    return lorentz_width, gauss_component

#%% Create macromolecule instances from JSON parameters
def create_MMs(tesla, params):
    # Calculate B0 is already done in load_mm_params; use positions directly
    names = params["names"]
    
    # Calculate linewidth components for both tissues
    lorentz_width_gm, gauss_component_gm = calc_lw_components(params["line_width_mms"], params["T2_gm"], 
                                                              params["line_width_water"], params["T2_water_gm"])
    
    lorentz_width_wm, gauss_component_wm = calc_lw_components(params["line_width_mms"], params["T2_wm"], 
                                                              params["line_width_water"], params["T2_water_wm"])
    
    # Normalize scales relative to the first peak
    scale_0 = params["scale"][0]
    scales = [round(s / scale_0, 4) for s in params["scale"]]

    # Calculate Bloch magnetization
    bloch_magnetization_gm = mm_bloch_simulation(params["T1_gm"], params["T2_gm"], plot=False)
    bloch_magnetization_wm = mm_bloch_simulation(params["T1_wm"], params["T2_wm"], plot=False)
    
    # Create Macromolecule instances
    instances = [
        Macromolecule(
            T1_pure_GM=params["T1_gm"][i],
            T1_pure_WM=params["T1_wm"][i],
            T2_pure_GM=params["T2_gm"][i],
            T2_pure_WM=params["T2_wm"][i],
            lorentz_width_gm=lorentz_width_gm[i],
            lorentz_width_wm=lorentz_width_wm[i],
            gauss_width_gm=gauss_component_gm[i],
            gauss_width_wm=gauss_component_wm[i],
            Bloch_magnetization_gm=bloch_magnetization_gm[i],
            Bloch_magnetization_wm=bloch_magnetization_wm[i],
            pos=params["position"][i],
            name=names[i],
            scale=scales[i],
            J_coupling_1=params["J1"][i],
            J_coupling_2=params["J2"][i]
        )
        for i in range(len(names))
    ]
    
    # Optionally, add each instance to globals if needed
    for mm in instances:
        globals()[mm.name] = mm
        
    return instances

#%% Simulation of MM spectra
def simulate_MMs(tesla, basis, TE, TR, json_file):
    """
    Simulate macromolecule spectra using the MM09 model.
    
    Args:
        tesla (float): Magnetic field strength in Tesla.
        basis (BasisSet): Basis set object containing spectral parameters.
        TE (float): Echo time in milliseconds.
        TR (float): Repetition time in milliseconds.
        json_file (str): Path to the JSON file containing macromolecule parameters.
    
    Returns:
        tuple: Simulated macromolecule spectra and areas under the MM092 peak.
    """
    params = load_mm_params(json_file, tesla)
    mm_peaks = create_MMs(tesla, params)

    tissues = ['WM', 'GM']
    mm_specs = torch.zeros([len(tissues), basis.n], dtype=torch.complex64)
    linspace = np.linspace(-basis.bw/2, basis.bw/2, basis.n)
    ppm_axis = basis.ppm
    time_axis = basis.t


    for tidx, tissue in enumerate(tissues):
        # Set tissue-specific parameters for each macromolecule
        for peak in mm_peaks:
            peak.set_T1(tissue)
            peak.set_T2_and_linewidths(tissue)
            peak.set_bloch_magnetization(tissue)
        
        # Vectorized initialization of MM lines
        mm_lines = np.array([voigt(peak.pos, peak.lorentz_width, peak.gauss_width, linspace) for peak in mm_peaks])
        
        # Apply all scaling factors to the MM lines
        for i, peak in enumerate(mm_peaks):
            # Scale with respect to MM092
            mm_lines[i, :] *= peak.scale
        
        # Correct for Bloch magnetization to obtain universal base MM spectrum
        mm_base = np.zeros_like(mm_lines)
        for i, peak in enumerate(mm_peaks):
            mm_base[i,:] = mm_lines[i, :] / peak.Bloch_magnetization

        # Sequence specific attenuation
        TE_sec = TE / 1000.0
        TR_sec = TR / 1000.0
       
        mm_spec = np.zeros_like(mm_base)
        for i, peak in enumerate(mm_peaks):
            if basis.localization == 'sLASER':
                att = Sequence.semiLASER(TE=TE_sec, TR=TR_sec, peak=peak)
            elif basis.localization == 'PRESS':
                att = Sequence.PRESS(TE=TE_sec, TR=TR_sec, peak=peak)
            else:
                raise ValueError(f"Unsupported localization: {basis.localization}")
            
            mm_spec[i, :] = att * mm_base[i, :]

        # Sum lines to form the real spectrum and compute Hilbert transform
        mm_spec = np.sum(mm_spec, axis=0)
        mm_spec = np.conjugate(hilbert(mm_spec))
        # Convert to a torch tensor
        mm_specs[tidx, :] = torch.tensor(mm_spec, dtype=torch.complex64)
   
    return mm_specs

