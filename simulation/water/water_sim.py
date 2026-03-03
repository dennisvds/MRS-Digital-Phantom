import numpy as np

#########################################
### Water Signal Simulation Functions ###
#########################################
# Based on the work of:
# Lin L, Považan M, Berrington A, Chen Z, Barker PB. 
# Water removal in MR spectroscopic imaging with L2 regularization. 
# Magnetic Resonance in Medicine. 2019;82(4):1278-1287. doi:10.1002/mrm.27824


def generate_sum_fid(time_axis, amplitudes, damping_hz, freq_ppm, phase_deg,
                     ref_freq_hz, decay_type='lorentzian', scaling=5.0):
    """
    Generate and sum multiple complex FIDs with either Lorentzian or Gaussian decay.

    Parameters
    ----------
    amplitudes : array_like, shape (M,)
        Signal amplitudes for each component.
    damping_hz : array_like, shape (M,)
        Linewidths (FWHM) in Hz for each component.
    freq_ppm : array_like, shape (M,)
        Precession frequencies in ppm for each component.
    phase_deg : array_like, shape (M,)
        Initial phases in degrees for each component.
    sw : float
        Spectral width (sampling rate) in Hz.
    n_points : int
        Number of time-domain points.
    ref_freq_hz : float
        Reference Larmor frequency in Hz (for ppm → Hz conversion).
    decay_type : {'lorentzian', 'gaussian'}
        Type of damping function (applied to all components).
    scaling : float
        Scaling factor (default is 5.0).

    Returns
    -------
    t : ndarray, shape (N,)
        Time vector (s).
    fid_sum : ndarray, shape (N,)
        Sum of all complex FID signals.
    """
    # Time axis
    t = time_axis

    # Convert inputs to arrays
    A = np.array(amplitudes)       # shape (M,)
    Γ = np.array(damping_hz)       # shape (M,)
    ppm = np.array(freq_ppm)       # shape (M,)
    φ = np.deg2rad(np.array(phase_deg))  # shape (M,)

    # Broadcast parameters: make shapes (M, 1)
    A = A[:, np.newaxis]
    Γ = Γ[:, np.newaxis]
    ppm = ppm[:, np.newaxis]
    φ = φ[:, np.newaxis]
    t_row = t[np.newaxis, :]       # shape (1, N)

    # Frequency conversion: ppm → Hz
    freq_hz = ppm * ref_freq_hz * 1e-6  # shape (M,1)

    # Damping
    if decay_type == 'lorentzian':
        alpha = np.pi * Γ
        decay = np.exp(-alpha * t_row)        # (M, N)
    elif decay_type == 'gaussian':
        beta = (np.pi * Γ)**2 / (4 * np.log(2))
        decay = np.exp(-beta * t_row**2)      # (M, N)
    else:
        raise ValueError("decay_type must be 'lorentzian' or 'gaussian'")

    # Complex precession for each component
    osc = np.exp(2j * np.pi * freq_hz * t_row + 1j * φ)  # (M, N)

    # Construct and sum FIDs
    fids = A * decay * osc         # (M, N)
    fid_sum = scaling * np.sum(fids, axis=0) # (N,)

    return fid_sum