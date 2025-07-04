from oct2py import Oct2Py
import os
import shutil
import json
import torch
import numpy as np
import hashlib
from scipy.io import loadmat, savemat

from gui.defaults import DEFAULT_METABOLITES

class MatlabRunner:
    def __init__(self, gui=None):
        """
        Initializes the MatlabRunner class.
        This class is responsible for running Octave simulations and generating basis sets.
        """
        self.octave = None
        self.gui = gui  # Optional GUI reference, can be used for progress updates or notifications
    
    def log(self, message: str, color: str = 'black'):
        """Helper function to log messages to GUI or print fallback."""
        if self.gui:
            self.gui.log_message(message, color)
        else:
            print(message)

    def start_octave(self):
        self.octave = Oct2Py()

        # Turn off warnings
        self.octave.eval("more off")
        self.octave.eval("echo off")
        self.octave.eval("warning('off', 'all')")

        # Add paths to the Octave session
        self.octave.addpath(self.octave.genpath('./simulation/FID-A'))
        self.octave.addpath(self.octave.genpath('./simulation/basissets/MRSCloud-main'))
        self.octave.addpath(self.octave.genpath('./simulation/lipids/SimnTG'))


    def stop_octave(self):
        """
        Stops the Octave session.
        """
        if hasattr(self, 'octave'):
            self.octave.exit()
            del self.octave
            self.log("Octave session stopped.")
        else:
            self.log("No Octave session to stop.", color='red')

    def generate_basis_set(self, settings):
        try:
            # Update simMRS.json file with the selected settings
            json_file_path = './simulation/basissets/MRSCloud-main/simMRS.json'
            with open(json_file_path, 'r') as file:
                json_data = json.load(file)

            # Modify json_data as needed
            json_data["userInput"]["vendor"] = [settings["vendor"]]
            json_data["userInput"]["localization"] = [settings["localization"]]
            json_data["userInput"]["TE"] = int(settings["TE"])
            json_data["userInput"]["metablist"] = DEFAULT_METABOLITES
            json_data["private"]["metab_default"] = []

            # Set save_dir
            json_data["DIRECTORY"]["save_dir"] = os.path.join(os.getcwd(), 'simulation', 'basissets', 'MRSCloud-main', 'basis_save')
            save_path = json_data["DIRECTORY"]["save_dir"]

            # Save the updated JSON data back to the file
            with open(json_file_path, 'w') as file:
                json.dump(json_data, file)

            # Run the Octave simulation
            self.log("Running Octave to generate basis set...")
            self.start_octave()
            self.octave.feval('run_simulations_cloud')

            # Move the generated basis set to the desired location
            basis_filename = f"LCModel_{settings['vendor']}_UnEdited_{settings['localization']}_TE{settings['TE']}.BASIS"
            source_path = os.path.join(save_path, basis_filename)
            destination_path = os.path.join('./data/basissets', basis_filename)

            # Ensure destination directory exists
            destination_dir = os.path.dirname(destination_path)
            if not os.path.exists(destination_dir):
                os.makedirs(destination_dir)

            # Check if the file exists before moving
            if os.path.exists(source_path):
                shutil.move(source_path, destination_path)
            else:
                print(f"Error: The file {source_path} does not exist.")

        except Exception as e:
            print(f"Error: {e}")
        finally:
            # Always exit Octave
            self.stop_octave()
    
    def generate_lipid_basis(self, basis, linewidth, components=True, cache_dir="./simulation/lipids/lipids_cache"):
        """
        Generates a lipid basis set using the Oct2Py interface.
        """
        # Call the lipid basis generation function in Octave

        # Parameters from basis set
        n_points = basis.n
        bw = 1 / basis.dwelltime
        B0 = basis.B0
        seq = basis.localization
        cf = basis.cf
        TE = basis.TE

        # Ensure cache directory exists
        os.makedirs(cache_dir, exist_ok=True)

        # Construct the path to the cache file
        key_string = f"{n_points}_{bw}_{B0}_{linewidth}_{TE}_{seq}_{cf}"
        hash_key = hashlib.md5(key_string.encode()).hexdigest()
        file_path = os.path.join(cache_dir, f"{hash_key}.mat")

        # Load from cache if it exists
        if os.path.isfile(file_path):
            self.log(f"Loading lipid basis set from cache: {file_path}")
            mat_data = loadmat(file_path)

            if components:
                fids = np.stack([mat_data[name] for name in mat_data if name.startswith("comp_")])
                fids = fids.reshape(len(fids), -1)
            else:
                fids = mat_data["fids"].reshape(1, -1)
        
        # If cache file does not exist, run Octave to generate the basis set
        else:
            self.log(f"Cache file not found. Generating lipid basis set: {file_path}")
            self.start_octave()
            out, out_components, names = self.octave.feval('create_lipids', n_points, bw, B0, linewidth, TE, seq, cf, nout=3)
            # Exit Octave
            self.stop_octave()
    
            # Save results
            if components:
                out_components_list = list(out_components[0])
                fids = []
                mat_dict = {}
                for i, struct in enumerate(out_components_list):
                    fid = struct['fids'][:, 0]
                    fids.append(fid)
                    mat_dict[f"comp_{i}"] = fid
                fids = np.array(fids)
            else:
                fid = out['fids'][:, 0]
                fids = np.array(fid).reshape(1, -1)
                mat_dict = {"fids": fids}
            
            # Save the generated basis set to cache
            savemat(file_path, mat_dict)
            self.log(f"Lipid basis set saved to cache: {file_path}")

      
        # Convert to complex64 (torch)
        lipid_fid = torch.tensor(np.array(fids), dtype=torch.complex64)

        return lipid_fid
