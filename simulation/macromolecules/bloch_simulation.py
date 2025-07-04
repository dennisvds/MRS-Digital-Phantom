import numpy as np
import matplotlib.pyplot as plt

def fSLR(B1, freqOffset, deltaT, rephasePoint):
    gamma_Hz = 42.576e6  # Hz/T
    gamma_rad = gamma_Hz * 2 * np.pi  # rad/(s*T)
    G = freqOffset / gamma_Hz  # Tesla

    a = np.ones((len(freqOffset), 2), dtype=complex)
    b = np.zeros((len(freqOffset), 2), dtype=complex)

    flagRephase = True

    for i in range(len(B1)):
        if flagRephase and i == rephasePoint:
            G = -G
            flagRephase = False

        B1eff = B1[i] * np.ones_like(freqOffset)
        phi = -gamma_rad * deltaT * np.sqrt(np.abs(B1eff)**2 + G**2)
        n0 = gamma_rad * deltaT / np.abs(phi)
        n = np.stack([
            np.real(B1eff) * n0,
            np.imag(B1eff) * n0,
            G * n0
        ], axis=1)

        av = np.cos(phi/2) - 1j * n[:,2] * np.sin(phi/2)
        bv = -1j * (n[:,0] + 1j * n[:,1]) * np.sin(phi/2)

        a[:,1] = av * a[:,0] - np.conj(bv) * b[:,0]
        b[:,1] = bv * a[:,0] + np.conj(av) * b[:,0]

        a[:,0], b[:,0] = a[:,1], b[:,1]

    return a[:,1], b[:,1]

def apply_rotation(a, b, M0):
    M1 = np.conj(a)**2 * M0[:,0] - b**2 * M0[:,1] + 2 * np.conj(a) * b * M0[:,2]
    M2 = -np.conj(b)**2 * M0[:,0] + a**2 * M0[:,1] + 2 * np.conj(b) * a * M0[:,2]
    M3 = -np.conj(a*b) * M0[:,0] - a*b * M0[:,1] + (np.abs(a)**2 - np.abs(b)**2) * M0[:,2]

    return np.stack([M1, M2, M3], axis=1).real

def apply_relaxation(M0, T1, T2, time):
    M = np.zeros_like(M0)
    M[:,0] = M0[:,0] * np.exp(-time / T2)
    M[:,1] = M0[:,1] * np.exp(-time / T2)
    M[:,2] = (M0[:,2] - 1) * np.exp(-time / T1) + 1
    return M



def mm_bloch_simulation(T1s, T2s, plot=False):
    # Convert T1 and T2 times to milliseconds
    T1s = np.array(T1s) * 1e3  # milliseconds
    T2s = np.array(T2s) * 1e3  # milliseconds

    # Time and pulse parameters
    time = 15e-3  # 15 ms
    nPoints = 512
    time_axis = np.linspace(0, time, nPoints) * 1e3  # in ms
    x = np.linspace(-np.pi, np.pi, nPoints)

    # AM and FM components of RF pulse
    b = 1.904
    AM = 1 / np.cosh(b * x)
    FM = 1017 * np.tanh(-x * 1.89)
    PM = 2 * np.pi * np.cumsum(FM) * time / len(FM)  # cumulative phase modulation

    # Frequency offsets
    freqOffset = np.arange(-2000, 2001, 10)  # Hz
    gammaHz = 42.576e6  # Hz/T

    # Initial magnetization M0 (along z)
    M0 = np.zeros((len(freqOffset), 3))
    M0[:, 2] = 1

    # Color map
    colors = plt.get_cmap('hsv', 10)

    # RF Pulse with some variation in amplitude
    ih = [0, 1, 2, -2, -1]
    B_pulses = [(1 + i * 0.15) * 15e-6 * AM * np.exp(1j * PM) for i in ih]  # list of complex B1 pulses

    # Use one RF pulse for both inversions (as in the last part of the script)
    Bampl = 24e-6  # T
    B = Bampl * AM * np.exp(1j * PM)

    # Compute alpha, beta
    deltaT = time / len(B)
    alpha, beta = fSLR(B, freqOffset, deltaT, rephasePoint=0)

    # DIR inversion times, DIR inversion times from Landheer et al
    TI1 = 920  # ms
    TI2 = 330  # ms

    mag_data = []

    for i in range(len(T1s)):
        T1 = T1s[i]
        T2 = T2s[i]

        # 1st inversion pulse
        M = apply_rotation(alpha, beta, M0)

        # Relaxation after first inversion
        M = apply_relaxation(M, T1, T2, TI1)

        # 2nd inversion pulse
        M = apply_rotation(alpha, beta, M)

        # Relaxation after second inversion
        M = apply_relaxation(M, T1, T2, TI2)

        # Store minimum Mz and plot
        mz = M[:, 2]
        mag_data.append(mz)
        

    # Final plot formatting
    if plot:
        plt.figure(figsize=(10, 6))
        for i in range(len(T1s)):
            mz = mag_data[i]
            plt.plot(freqOffset, mz, color=colors(i), label=f'M{i+1}', linewidth=1)
        plt.xlabel('Frequency (Hz)', fontsize=15, fontweight='bold')
        plt.ylabel(r'$M_z/M_0$', fontsize=15, fontweight='bold')
        plt.title('Bloch Simulation of MM DIR Response', fontsize=14)
        plt.gca().invert_xaxis()
        plt.legend(fontsize=9, loc='lower right')
        plt.grid(True)
        plt.axis([-2000, 2000, -1, 1])
        plt.tight_layout()
        plt.show()

    mag_data = np.array(mag_data)
    bloch_magnetization = np.min(mag_data, axis=1)

    return bloch_magnetization


if __name__ == "__main__":
    # Example usage GM
    T1s = np.array([290, 309, 309, 225, 247, 263, 400, 400, 400, 400])* 1e-3  # seconds
    T2s = np.array([27, 39, 18, 17, 14, 20, 26, 21, 18, 21])* 1e-3  # seconds

    # Example usage WM
    # T1s = np.array([290, 309, 309, 225, 247, 263, 400, 400, 400, 400])* 1e-3  # seconds
    # T2s = np.array([27, 39, 18, 17, 14, 20, 26, 21, 18, 21])* 1e-3  # seconds

    bloch_magnetization = mm_bloch_simulation(T1s, T2s, plot=True)
    print("Bloch magnetization values:", bloch_magnetization)

    # GM: [0.38589514 0.34759162 0.34759162 0.54634556 0.4869087  0.44696159 0.211404   0.211404   0.211404   0.211404  ]
    # WM: [0.39878293 0.40758989 0.43036719 0.47160023 0.32892172 0.34000284 0.23153754 0.12138586 0.23659643 0.17726085]

    # Current values (No difference in GM and WM):
    # GM: [0.38589514 0.34759162 0.34759162 0.54634556 0.4869087  0.44696159 0.211404   0.211404   0.211404   0.211404  ]
    # WM: [0.38589514 0.34759162 0.34759162 0.54634556 0.4869087  0.44696159 0.211404   0.211404   0.211404   0.211404  ]