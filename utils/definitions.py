####################################################################################################
#                                       Metabolite names                                           #
####################################################################################################
METABS = [
    'AcAc',   # Acetoacetate
    'Ace',    # Acetate
    'Ac0',    # Acetone
    'Ala',    # Alanine
    'Asc',    # Ascorbate
    'Asp',    # Aspartate
    'Bet',    # Betaine
    'bHB',    # beta-hydroxybutyrate
    'bHG',    # 2-hydroxyglutarate
    'Car',    # Carnitine
    'Cit',    # Citrate
    'Cr',     # Creatine
    'Cys',    # Cysteic acid
    'Cystat', # Cystat
    'CrCH2',  # negative CrCH2 correction signal
    'EA',     # Ethanolamine
    'EtOH',   # Ethanol
    'fCho',   # free choline
    'Fuc',    # Fucose
    'GABA',   # GABA
    'Gcn',    # Glucone
    'Gcr',    # Glucoronic acid
    'GPC',    # Glycerophosphocholine
    'GSH',    # Glutathione (reduced)
    'Glc',    # Glucose
    'Gln',    # Glutamine
    'Glu',    # Glutamate
    'Gly',    # Glycine
    'Gua',    # Guanidinoacetate
    'H2O',    # H2O
    'HCar',   # Homocarnosine
    'ILc',    # Isoleucine
    'mI',     # myo-inositol
    'Lac',    # Lactate
    'Leu',    # Leucine
    'Lys',    # Lysine
    'NAA',    # N-Acetylaspartate
    'NAAG',   # N-Acetylaspartylglutamate
    'PCh',    # phosphorylcholine
    'PCr',    # Phosphocreatine
    'PE',     # Phosphoethanolamine
    'Pgc',    # Propyleneglycol
    'Phenyl', # Phenylalanine
    'Pyr',    # Pyruvate
    'sI',     # scyllo-inositol
    'Ser',    # Serine
    'Suc',    # Succinate
    'Tau',    # Taurine
    'Thr',    # Threonine
    'Tyros',  # Tyrosine
    'Val',    # Valine
    'NAA_Ace',# NAA acetyl
    'NAA_Asp',# NAA aspartyl

    # Macromolecules (MM + ppm value)
    'MM09',   # MM09
    'MM12',   # MM12
    'MM14',   # MM14
    'MM17',   # MM17
    'MM20',   # MM20
    'MM22',   # MM22
    'MM26',   # MM26
    'MM27',   # MM27
    'MM30',   # MM30
    'MM32',   # MM32
    'MM36',   # MM36
    'MM37',   # MM37
    'MM38',   # MM38
    'MM42',   # MM42

    # Reference at 0 ppm
    'Ref0ppm',
]

####################################################################################################
#                             Metabolite translation dictionary                                    #
####################################################################################################
METABS_TRANSLATION = {
    # Creatine
    'Cre': 'Cr',
    'Cre_CH2': 'Cr',
    'Cre_CH3': 'Cr',
    # Choline
    'Cho': 'fCho',
    # Myo-inositol
    'Myo': 'mI',
    'Ins': 'mI',
    # Scyllo-inositol
    'Scy': 'sI',
    'Scyllo': 'sI',
    # Glucose
    'Glc ': 'Glc',
    # 'Glc_B': 'Glc', # Glucose B (not in the list)
    # Glycerophosphocholine
    'GPC_Gly': 'GPC',
    'GPC_Mult': 'GPC',
    # N-Acetylaspartylglutamate
    'NAAG_CH3': 'NAAG',
    'NAAG_Asp': 'NAAG',
    'NAAG_Glu': 'NAAG',
    # N-Acetylaspartate
    'NAA_CH3': 'NAA',
    'NAA_Asp': 'NAA',
    # Phosphocholine
    'PCho': 'PCh',
}

####################################################################################################
#                                       Label names BigBrain-MR                                    #
####################################################################################################
LABELS = {
    0: "Background",
    1: "White Matter",
    2: "Gray Matter",
    3: "Cerebrospinal Fluid",
    4: "Cerebellum White Matter",
    5: "Cerebellum Gray Matter",
    6: "Thalamus",
    7: "Caudate",
    8: "Putamen",
    9: "External Pallidum",
    10: "Internal Pallidum",
    11: "Basal Forebrain",
    12: "Accumbens",
    13: "Brainstem",
    14: "Hippocampus",
    15: "Ventral Diencephalon",
    16: "Amygdala",
    17: "Subthalamic Nuclei",
    18: "Red Nuclei",
    19: "Substantia Nigra",
    20: "Pineal Gland",
}

####################################################################################################
#                                     Gyro-Magnetic Ratios                                         #
####################################################################################################
# In [MHz/T]
GYRO_MAG_RATIO = {
    '1H': 42.576,
    '2H': 6.536,
    '13C': 10.7084,
    '31P': 17.235}

####################################################################################################
#                                           CSF data                                               #
####################################################################################################
# 'metabolite': [conc_mean, conc_std, T1, T2]
# Literature:
# [1]   Minati, L., Aquino, D., Bruzzone, M. G., & Erbetta, A. (2010). 
#       Quantitation of normal metabolite concentrations in six brain regions by in-vivo 1H-MR spectroscopy. 
#       Journal of Medical Physics / Association of Medical Physicists of India, 35(3), 154–163. https://doi.org/10.4103/0971-6203.62128

CSF_DATA = {
    'H2O': [55.5, 0, 0, 160], # Based on [1] (T2 based on 1.5T, assuming same for 3T)
    'Lac': [0.1, 0, 0, 110],    # Manual input, same T2 as GM in metabolite dataframe
}

#####################################################################################################
#                                       Water Concentrations                                        #
#####################################################################################################
# Based on literature values:
# Gasparovic, C., Neeb, H., Feis, D.L., Damaraju, E., Chen, H., Doty, M.J., South, D.M., Mullins, P.G., Bockholt, H.J. and Shah, N.J. (2009), 
# Quantitative spectroscopic imaging with in situ measurements of tissue water T1, T2, and density. 
# Magn. Reson. Med., 62: 583-590. https://doi.org/10.1002/mrm.22060

H2O_CONCENTRATIONS = {
    'Background': 0.0,
    'WM': 0.703 * 55.5,
    'GM': 0.851 * 55.5,
    'CSF': 55.5,
    'Fat': 0.0,
    'Skull': 0.0,
}