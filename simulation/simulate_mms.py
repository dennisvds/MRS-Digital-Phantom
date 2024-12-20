####################################################################################################
#                                      simulate_mms.py                                             #
####################################################################################################
#                                                                                                  #
# Authors: Saipavitra Murali-Manohar & Andrew Martin Wright                                        #
#          adapted by: D.M.J. van de Sande (d.m.j.v.d.sande@tue.nl)                                #
#                                                                                                  #
# Created: 03/02/23                                                                                #
#                                                                                                  #
# Purpose: Functions and classes used to create MM spectra. With the goal being to use the spectra #
#           as MM baselines in metabolite spectra.                                                 #
#                                                                                                  #                
# This code is part of the work:                                                                   #                
# Wright AM, Murali-Manohar S, Borbath T, Avdievich NI, Henning A.                                 #
# Relaxation-corrected macromolecular model enables determination of 1H longitudinal T1-relaxation #
# times and concentrations of human brain metabolites at 9.4T.                                     #      
# Magnetic Resonance in Medicine. 2022;87(1):33-49. doi:10.1002/mrm.28958                          #
#                                                                                                  #  
####################################################################################################

import numpy as np


#%% Macromolecule class for defining macromolecules
# might be useful in including a metabolite class
class Macromolecule():
    """
    Macromolecule characteristics to be used in attenuation calculations.
    Characteristics include T1, T2, position along the index axis, the 
    name of the peak, and the scale of the peak from a DIR_semiLASER sequence 
    relative to MM09.
    """
    def __init__(self, T1_pure_GM, T1_pure_WM, T2, T2_star, lorentz_width, gauss_width, Bloch_magnetization, pos, name, scale, J_coupling_1, J_coupling_2):
        self.T1_GM = T1_pure_GM
        self.T1_WM = T1_pure_WM
        self.T2 = T2
        self.T2_star = T2_star
        self.pos = pos
        self.name = name
        self.scale = scale #Relative to MM09 from a DIR-semiLASER sequence
        self.lorentz_width = lorentz_width
        self.Bloch_magnetization = Bloch_magnetization
        self.gauss_width = gauss_width
        self.J1 = J_coupling_1
        self.J2 = J_coupling_2
        self.T1 = np.empty(1)
    
    def add_voxel_T1(self, T1):
#        self.T1.append(T1)
        self.T1 = T1

#%% creation of Voigt lines
def voigt(pos, g, s, x):
    """
    creates voigt lines for the MM model
    x = 1:4096 number of points
    pos = offset for the center of the MM peak in points
    g = FWHM of Lorentzian
    s = FWHM of Gaussian
    These voigt shapes are further normalized to all have amplitude = 1.
    """
    pos = pos/2.0
    s = s/np.sqrt(8*np.log(2))
    gauss = (1/(s*np.sqrt(2*np.pi))) * np.exp(-(((x-pos))**2)/(2*s**2))
#    gauss = (1/(s*np.sqrt(2*np.pi))) * np.exp(-(((pos-x))**2)/(2*s**2))
    
#    g = 1/(np.pi * g)
    #a = 1 + (((x-pos/2))/(g))**2
    #lorentz = (1/(np.pi*g)) * (1/a)
    a = g/2.0
    lorentz = (1/np.pi) * ((a)/((x-pos)**2 + a**2))
    line = np.convolve(gauss, lorentz, mode = 'same')
    norm_line = line/np.max(line)
    
    return norm_line

#%% defining sequences to be used for attenuation
class sequence():
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

######################
### Phantom inputs ###
######################
# Values are based on:
#
# [1]   Landheer K, Gajdošík M, Treacy M, Juchem C. Concentration and effective T2 relaxation times of macromolecules at 3T. 
#       Magnetic Resonance in Medicine. 2020;84(5):2327-2337. doi:10.1002/mrm.28282.
#
# [2]   Hoefemann M, Bolliger CS, Chong DGQ, van der Veen JW, Kreis R. Parameterization of metabolite and macromolecule 
#       contributions in interrelated MR spectra of human brain using multidimensional modeling. NMR in Biomedicine. 2020;33(9):e4328. doi:10.1002/nbm.4328
#       as reported in:
#       Cudalbu C, Behar KL, Bhattacharyya PK, et al. Contribution of macromolecules to brain 1H MR spectra: Experts’ consensus recommendations. 
#       NMR in Biomedicine. 2021;34(5):e4393. doi:10.1002/nbm.4393
# 

def create_MMs(tesla):
    B0    = 42.576 * tesla #Hz
    name = ["MM092", "MM121", "MM139", "MM167", "MM204", "MM226", "MM270", "MM299", "MM321", "MM375"] # from [1], ppm naming convention
    # T1-vales
    T1_gm = np.array([290, 309, 309, 225, 247, 263, 400, 400, 400, 400]) / 1000     # from [2], in seconds
    T1_wm = np.zeros_like(T1_gm)                                                    # not defined
    # T2-values
    T2_gm = np.array([26.6, 38.7, 17.9, 16.6, 14.3, 19.9, 25.5, 21.0, 17.7, 21.0]) / 1000 # from [1], in seconds
    T2_wm = np.zeros_like(T2_gm)                                                          # not defined
    T2_gm_water = 79.9 / 1000                                                                  
    # T2*-values
    T2_star_gm = np.zeros_like(T2_gm) # not defined
    T2_star_wm = np.zeros_like(T2_wm) # not defined
    # Other values
    lorentz_width_water = 1/(np.pi * T2_gm_water)
    shim_component = np.ones(len(name)) * 7.9                                                                                           # avg water linewidth from [1] in Hz
    line_width = np.array([(25.2),(21.1),(28.5),(38.2),(37.1),(32.4),(35.3),(28.0),(17.2),(57.0)])                                      # LW from [1] in Hz
    Bloch_magnetization = [0.6455, 0.6724, 0.6442, 0.5251, 0.4246, 0.5317, 0.5545, 0.6207, 0.3627, 0.5641, 0.1297, 0.4035, 0.6675]      # from original input
    position = np.array([0.92*B0, 1.21*B0, 1.39*B0, 1.67*B0, 2.04*B0, 2.26*B0, 2.70*B0, 2.99*B0, 3.21*B0, 3.75*B0]) - (4.65*B0)         # ppm values from [1]
    scale = [21.1, 8.2, 20.4, 48.0, 78.4, 50.4, 19.0, 29.1, 10.0, 40.1]                                                                 # concentrations from Pavi
    J1 = [7.3, 6.9, 7.3, 0, 0, 0, -17.5, 0, 0, 0]                                                                                       # from original input
    J2 = [0, 0, 0, 0, 0, 0, 6, 0, 0, 0]                                                                                                 # from original input

    ###################
    ### Computation ###
    ###################
    lorentz_width = 1/ (np.pi * T2_gm) #in Hz
    overlap_component = line_width - lorentz_width - shim_component 
    gauss_component = shim_component + overlap_component
    gauss_component = np.clip(gauss_component, a_min=0.1, a_max=None) #in Hz

    scale_0 = scale[0]
    for idx, values in enumerate(scale):
        scale[idx] = round(values / scale_0, 4)


    instances = [Macromolecule(T1_gm[idx], T1_wm[idx], T2_gm[idx], T2_star_gm[idx], lorentz_width[idx], gauss_component[idx], Bloch_magnetization[idx], position[idx], name[idx], scale[idx], J1[idx], J2[idx]) for idx, mm in enumerate(name)]

    for idx, mm in enumerate(name):
        globals()[f"{mm}"] = instances[idx]

    MM = instances

    return MM

def simulate_MMs(tesla, basis, labels, scaling):
    MM = create_MMs(tesla)

    mm_peaks = np.array(MM)

     # Set T1 to GM T1
    for i in range(len(mm_peaks)):
        mm_peaks[i].T1 = mm_peaks[i].T1_GM

    # Initialize MM lines
    mm_lines = np.zeros([len(mm_peaks), basis.n])

    # Scaling
    linspace = np.linspace(-basis.bw/2, basis.bw/2, basis.n)
    for i in range(len(mm_peaks)):
        mm_lines[i, :] = voigt(mm_peaks[i].pos, mm_peaks[i].lorentz_width, mm_peaks[i].gauss_width, linspace) * mm_peaks[i].scale
    
    # Magnetization
    for i in range(len(mm_peaks)):
        mm_lines[i, :] = mm_lines[i, :] / mm_peaks[i].Bloch_magnetization

    # Attenuation
    attenuation = np.zeros([len(mm_peaks), 1])
    for i in range(len(mm_peaks)):
        attenuation[i] = sequence.FID(TE=0.025, TR=1.0, peak=mm_peaks[i])
    
    # Resultant spectrum
    resultant_spectrum = np.dot(mm_lines.T, attenuation).T

    # Reshape to the desired output shape
    resultant_spectrum = np.broadcast_to(resultant_spectrum, (len(labels), basis.n))

    # Create a writable copy of resultant_spectrum
    resultant_spectrum = np.copy(resultant_spectrum)

    # Remove MM baseline from all labels other than GM (1) and WM (2)
    resultant_spectrum[(labels != 1) & (labels != 2)] = 0

    return resultant_spectrum * scaling  

# Notes:
    # Old scaling:
    # MM spectrum is divided by 12 to match the amplitude of the metabolite signals. 
    # 12 is based on the number of protons in the MM 0.92 peak (Isoleucine, L-Leucine, L-Valine), since scale factors were based on this peak.
    # Example Isoleucine: https://hmdb.ca/spectra/nmr_one_d/1136
    # Scaling is now done with a scaling factor in the overall signal model to control the MM amplitudes.