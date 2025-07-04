import os

# paths
DEFAULT_METAB_DF_PATH = "./data/metabolites/metab_df.csv"
DEFAULT_MM_JSON_PATH = "./data/macromolecules/mm_params.json"
DEFAULT_SKELETON_DIR = "./data/skeletons/"
DEFAULT_BASIS_SET_DIR = "./data/basissets/"


# metabolite list used for simulation (MRSCloud), based on FID-A spin systems available
# Removed from FID-A defaults:
# "Lys" - Lysine -- simulation time too long
# "HCar" - H-Creatine -- spin system not available in FID-A spinSystems.mat
# "2HG" - 2-Hydroxyglutarate -- Error in compiling .basis file
# "NAA_Aspartate_only" - NAA Aspartate -- Not needed for simulation
# "NAA_Acetyl_only" - NAA Acetyl -- Not needed for simulation

DEFAULT_METABOLITES = ['H2O', 'Tau', 'Scyllo', 'Gln', 'Gly', 'bHB', 'PCh', 'bHG', 'Lac', 'NAA', 'GPC', 
                       'GSH', 'Cit', 'Asc', 'GABA', 'Thr', 'EtOH', 'NAAG', 'Asp', 'Cr', 'Phenyl', 'Ser', 
                       'Ref0ppm', 'EA', 'Ala', 'PCr', 'Glc', 'Ins', 'PE', 'Glu', 'Val', 'Tyros']

# skeleton options
DEFAULT_SKELETON = "MRiLab"

# simulation defaults
DEFAULT_ECHO_TIME_MS = 35
DEFAULT_SPECTRAL_POINTS = 2048
DEFAULT_REPETITION_TIME_MS = 2000
DEFAULT_NOISE_LEVEL = 5.0
DEFAULT_MM_LEVEL = 25.0
DEFAULT_BW = 2000
DEFAULT_VOXEL_SIZE_MM = 10.00

# lipid defaults
DEFAULT_LIPID_AMP_FACTOR = 40.0
DEFAULT_LIPID_SIGMA = 5.0
DEFAULT_LIPID_PHASE_MIN = -180 # degrees
DEFAULT_LIPID_PHASE_MAX = 180 # degrees
DEFAULT_LIPID_LW_MIN = 30.0 # Hz
DEFAULT_LIPID_LW_MAX = 50.0 # Hz

# Water defaults
DEFAULT_WATER_AMP_FACTOR = 5.0
DEFAULT_WATER_PHASE_MIN = -180 # degrees
DEFAULT_WATER_PHASE_MAX = 180 # degrees
DEFAULT_WATER_DAMPING_MIN = 0 # Hz
DEFAULT_WATER_DAMPING_MAX = 50 # Hz

# Shim imperfection defaults
DEFAULT_SHIM_AMPLITUDE_HZ = 3.0 # Hz
DEFAULT_SHIM_CORR_LENGTH = 1.0 # mm
DEFAULT_SHIM_BOUNDARY_AMP_FACTOR = 1.5
DEFAULT_SHIM_BOUNDARY_SMOOTHING = 2.0 # mm

# Defaults for batch simulations
DEFAULT_path2basis = os.path.join(DEFAULT_BASIS_SET_DIR, [f for f in os.listdir(DEFAULT_BASIS_SET_DIR) if f.endswith('.BASIS')][0])
DEFAULT_SIM_PARAMS = {
            "spectral_points": DEFAULT_SPECTRAL_POINTS,
            "TR": DEFAULT_REPETITION_TIME_MS,
            "bandwidth": DEFAULT_BW,
            "noise_level": DEFAULT_NOISE_LEVEL,
            "mm_level": DEFAULT_MM_LEVEL,
            "lipid_amp_factor": DEFAULT_LIPID_AMP_FACTOR,
            "lipid_sigma": DEFAULT_LIPID_SIGMA,
            "shim_amplitude_hz": DEFAULT_SHIM_AMPLITUDE_HZ,
            "shim_corr_length": DEFAULT_SHIM_CORR_LENGTH,
            "shim_boundary_amp_factor": DEFAULT_SHIM_BOUNDARY_AMP_FACTOR,
            "shim_boundary_smoothing": DEFAULT_SHIM_BOUNDARY_SMOOTHING,
            "water_amp_factor": DEFAULT_WATER_AMP_FACTOR,
            "TE": DEFAULT_ECHO_TIME_MS,
            "vendor": "Philips",
            "localization": "PRESS",
        }