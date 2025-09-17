# MRS Digital Brain Phantom Framework

## Overview
This repository contains the code for the **MRS Digital Brain Phantom Framework**, a tool for simulating Magnetic Resonance Spectroscopy (MRS) data. The framework is designed to be comprehensive and flexible, allowing researchers to generate MRS data efficiently for their own applications.

<figure style="text-align: center;">
    <img src="media/framework.png" alt="Sample Figure" width="800" />
    <figcaption>Figure 1: Structure of the 3D MRS digital brain phantom framework. The framework is divided into three stages: Skeleton, MRS Phantom, and Simulation. Each stage allows for user-defined inputs, making the framework highly modular and customizable for various applications. This example uses the ‘BigBrain-MR’ skeleton source.</figcaption>
</figure>

## Getting Started

### 1. Downloading Data
To use the **BigBrain-MR phantom**, you can download the required data [**here**](https://zenodo.org/records/7432527).  
The following files are needed for this project:  
- `BigBrainMR_Labels_400um.nii.gz`  
- `BigBrainMR_T1map_400um.nii.gz`  
- `BigBrainMR_T2smodel_R2smap_400um.nii.gz`  
- `BigBrainMR_Tweighted_400um.nii.gz`  

All files should be placed in the folder: `data/skeletons/BigBrainMR/`.  

All other required data are included in this repository. If you would like to use your own skeletons, some adjustments to the code may be required.  


### 2. Installation
You can install and run the simulator in two ways:
- **Option A:** Pixi (for Mac/Linux)
- **Option B:** Conda (cross-platform)

### Option A: Install & Run with Pixi

On macOS and Linux you can use the [Pixi](https://pixi.sh/latest/) package manager to run the simulator.

First [install Pixi](https://pixi.sh/latest/installation/) by running `curl -fsSL https://pixi.sh/install.sh | sh`

Then, clone this repository:
```bash
git clone https://github.com/dennisvds/MRS-Digital-Phantom.git
cd MRS-Digital-Phantom
```

Then run the simulator by calling `pixi run simulator` in the cloned repository.

### Option B: Install & Run with Conda
#### Prerequisites

Ensure the following dependencies are installed:

- **[Python ≥ 3.10](https://www.python.org/)**  
- **[Conda](https://docs.conda.io/en/latest/)**  
- **[Jupyter](https://jupyter.org/install)** (JupyterLab or Jupyter Notebook, to run `.ipynb` files)

- **GNU Octave ≥ 4.0**  
  Required for simulating lipid signals and basis sets.

  **Installation instructions for GNU:**

  - **Windows**  
    Download the latest *MinGW* version from the [GNU Octave website](https://www.gnu.org/software/octave/download.html).  
    During installation, make sure to add Octave to your system's `PATH`. You can do this during setup or by manually adding the path to the `bin` directory to your environment variables.

  - **macOS**  
    Install using Homebrew:

    ```bash
    brew install octave
    ```

  - **Linux**  
    Use your distribution’s package manager. For example, on Ubuntu:

    ```bash
    sudo apt-get update
    sudo apt-get install octave
    ```

  **To verify installation**, run:

  ```bash
  octave --version

Continue the installation by doing the following steps:

#### 1. Clone the repository

```bash
git clone https://github.com/dennisvds/MRS-Digital-Phantom.git
cd MRS-Digital-Phantom
```

#### 2. Create and activate a new virtual environment 
```bash
conda create -n mrs-phantom python=3.10
conda activate mrs-phantom
```

#### 3. Install dependencies
```bash
conda install -c conda-forge -c https://fs/.fmrib.ox.ac.uk/fsldownloads/fslconda/public/
fsl_mrs octave ipython jupyter_core pyqtgraph pytorch oct2py torchio
```

#### 4. Run the simulator
With everything installed and the data in the correct folder structure, you can run the MRS digital phantom:
```bash
python main.py
```

## Directory Structure
The following directory structure is used in this project:
```
MRS_Digital_Phantom/
├── data/                           # Main folder to store the data
│   ├── basissets/                  # Folder to store basissets
│   ├── invivo/                     # In-vivo data for comparison analysis
│   ├── macromolecules/             # Data for metabolite
│   ├── metabolites/                # Data from the MRiLab phantom
│   ├── skeletons/                  # Folder to store all MRS Phantom data
├── gui/                            # Scripts needed to generate and run the GUI elements
├── installation/                   # Scripts for installing the necessary packages
├── media/                          # Folder containing images for illustration
├── notebooks/                      # Folder containing notebooks used for results in the paper
│   ├── analysis.ipynb              # Notebook used to generate paper results for in-vivo vs simulated data comparison
│   ├── generate_data.ipynb         # Notebook used to generate batches of simulated data using configuration files
│   ├── plot_variations.ipynb       # Notebook used to generate paper results for variations in simulated data
├── outputs/                        # Folder containing simulation outputs used in the publication
├── simulation/                     # Folder with simulation scripts and signal models
│   ├── basissets/                  # Code for generating basissets (using MRSCloud)
│   ├── FID-A/                      # FID-A code needed for basisset and lipid generation
│   ├── lipids/                     # Code for generating lipid signal (based on SimnTG)
│   ├── macromolecules/             # Code for generating macromolecules signals
│   ├── water                       # Code for generating residual water signals
├── utils/                          # Utility scripts for plotting, loading, and setting definitions
├── config.json                     # Example of a configuration file in .JSON format
├── main.py                         # Main script that runs the GUI
├── environment.yml                 # Conda environment file
├── requirements.txt                # Python dependencies
└── README.md                       # Project documentation (this file)
```

## Usage
When the GUI is opened, you will see the following screen:

<figure style="text-align: center;">
    <img src="media/GUI.png" alt="Sample Figure" width="800" />
    <figcaption> Figure 2: Screenshot of the graphical user interface (GUI) of the digital MRS phantom. The left panel (green box) displays settings for the skeleton, metabolite dataframe, and basis set. The top panel (purple box) shows the three orthogonal brain views used for voxel placement. The bottom-middle section (blue box) contains the simulation settings panel and message box. The right panel (red box) visualizes the simulated spectrum, including its individual signal components. </figcaption>
</figure>

> #### Stability Notes
> While the GUI is fully functional and has been tested across platforms, users may occasionally encounter instability (e.g., unexpected crashes or segmentation faults), especially when using certain Python/Qt configurations or interactive environments. We recommend launching the GUI from a clean virtual environment and avoiding conflicting toolkits or background processes. > If issues persist, restarting the environment or switching between `pip` and `conda` installations of PyQt5 may help.
>
> In addition, the **BigBrain-MR phantom** requires substantial memory due to its high resolution. Simulations with moderate/large voxel sizes (> 2.5 cm) can easily exceed the available RAM on standard laptops or desktops and may cause crashes. We recommend either using smaller voxel sizes for testing or running BigBrain simulations on machines with ample memory.  
 

### Notebooks
There are a couple of Jupyter Notebooks available in this repository that will go over the batch-wise generation of spectra and analysis methods. These notebooks were used to create the figures in the corresponding paper. [View the Notebooks](https://github.com/dennisvds/MRS_Digital_Phantom/blob/main/Demo_MRS_Phantom.ipynb).

## Contact, Feedback, Suggestions
This MRS Digital Brain Phantom is a never-ending project! We are happy to receive your questions, feedback, suggestions, or critique by sending an email to:
**Dennis van de Sande** (d.m.j.v.d.sande@tue.nl)

## Citation
Should you publish any material that made use of this MRS Phantom Framework, please cite the following paper:

*Preprint (currently under review)*:

[D.M.J. van de Sande, A.T. Gudmundson, S. Murali-Manohar, C.W. Davies-Jenkins, D. Simicic, G. Simegn, İ. Özdemir, S. Amirrajab, J.P. Merkofer, H.J. Zöllner, G. Oeltzschner, R.A.E. Edden A Digital Phantom for MR Spectroscopy Data Simulation.
](https://arxiv.org/abs/2412.15869)


## License
This project is licensed under the [MIT License](https://github.com/dennisvds/MRS_Digital_Phantom/blob/main/LICENSE.md).

## Acknowledgments
This work has been (partially) funded by:
- European ITEA4 program (project 20209)