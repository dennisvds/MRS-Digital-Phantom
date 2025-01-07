# MRS Digital Brain Phantom Framework

## Overview
This repository contains the code for the **3D MRS Digital Brain Phantom Framework**, a tool for simulating Magnetic Resonance Spectroscopy (MRS) and Spectroscopic Imaging (MRSI) data. The framework is designed to be comprehensive and flexible, allowing researchers to generate MRS data efficiently for their own applications.

<figure style="text-align: center;">
    <img src="media/framework.png" alt="Sample Figure" width="800" />
</figure>

## Getting Started

### Prerequisites
Make sure you have the following installed:

- [Python](https://www.python.org/) (>= 3.9)
- [conda](https://docs.conda.io/en/latest/) (optional, but recommended)
- [Jupyter](https://jupyter.org/install) (JupyterLab or Jupyter Notebook to run the .ipynb files)

### Installation

#### Option 1: Using `conda` Environment
1. Clone the repository:
   ```bash
   git clone https://github.com/dennisvds/MRS_Digital_Phantom.git
   cd MRS_Digital_Phantom
   ```
2. Create an environment using conda:
   ```bash
   conda create --name brainphantom python=3.9
   ```
   *Note:* Python version should be **3.9** or higher.

3. Activate the environment:
   ```bash
   conda activate brainphantom
   ```
4. Install packages using script:
   ```bash
   ./installation/install_conda.sh
   ```

#### Option 2: Using Python Virtual Environment (`venv`)
1. Clone the repository:
   ```bash
   git clone https://github.com/dennisvds/MRS_Digital_Phantom.git
   cd MRS_Digital_Phantom
   ```
2. Create a virtual environment:
   ```bash
   python -m venv brainphantom
   ```
3. Activate virtual environment:
   on macOS/Linux:
   ```bash
   source brainphantom/bin/activate
   ```
   on Windows:
   ```bash
   brainphantom\Scripts\activate
   ```
4. Install packages using script:
   ```bash
   ./installation/install_pip.sh
   ```
   *Note:* FSL-MRS does not have a PIP-installation. Refer to the [FSL-MRS Documentation](https://open.win.ox.ac.uk/pages/fsl/fsl_mrs/install.html) to install this package.

#### Packages
If you want install the packages yourself, here is a list of the packages that are installed:
- ipykernel
- fsl_mrs
- torchio
- ipywidget
- openxyl
- ipympl

Development has been done using Python 3.12

### Data
Data for the **BigBrain-MR phantom** can be downloaded [**here**](https://zenodo.org/records/7432527). All other data is available in this respository. If you want to add your own skeletons and/or metabolite dataframe, you can add or change the files in the data folder. 

### Example Usage
There is a Demonstration Jupyter Notebook in this repository that will go over all steps of generating the MRS Phantom and on how to generate spectral data using it. [View the Notebook](https://github.com/dennisvds/MRS_Digital_Phantom/blob/main/Demo_MRS_Phantom.ipynb).


## Directory Structure
The following directory structure is used in this project:
```
MRS_Digital_Phantom/
├── data/                           # Main folder to store the data
│   ├── Bassissets/                 # Folder to store basissets
│   ├── BigBrainMR/                 # Data from the BigBrainMR phantom
│   ├── metabolites/                # Data for metabolite
│   ├── MRiLab/                     # Data from the MRiLab phantom
│   ├── phantom/                    # Folder to store all MRS Phantom data
├── installation/                   # Scripts for installing the necessary packages
├── loading/                        # Scripts for loading files
├── media/                          # Folder containing images for illustration
├── preprocessing/                  # Preprocessing scripts
├── simulation/                     # Simulation scripts and signal models 
├── utils/                          # Utility scripts
├── computation_test.ipynb          # Notebook used to measure computation times
├── create_examples.ipynb           # Notebook that shows different variations of spectral simulations
├── create_figures.ipynb            # Notebook that is used to generate images for the paper
├── Demo_MRS_Phantom.ipynb          # Main Notebook that showcases the MRS Phantom usage
├── visualize_mrsi.ipynb            # Notebook to generate MRSI plots (used for paper)
├── environment.yml                 # Conda environment file
├── requirements.txt                # Python dependencies
└── README.md                       # Project documentation (this file)
```


## Contact, Feedback, Suggestions
This MRS Digital Brain Phantom is a never-ending project! We are happy to receive your questions, feedback, suggestions, or critique by sending an email to:
**Dennis van de Sande** (d.m.j.v.d.sande@tue.nl)

## Citation
Should you publish any material that made use of this MRS Phantom Framework, please cite the following paper:

*Preprint (currently under review)*:

[D.M.J. van de Sande, A.T. Gudmundson, S. Murali-Manohar, C.W. Davies-Jenkins, D. Simicic, G. Simegn, İ. Özdemir, S. Amirrajab, J.P. Merkofer, H.J. Zöllner, G. Oeltzschner, R.A.E. Edden A Digital Phantom for 3D MR Spectroscopy Data Simulation.
](https://arxiv.org/abs/2412.15869)


## License
This project is licensed under the [MIT License](https://github.com/dennisvds/MRS_Digital_Phantom/blob/main/LICENSE.md).

## Acknowledgments
This work has been (partially) funded by:
- European ITEA4 program (project 20209)
