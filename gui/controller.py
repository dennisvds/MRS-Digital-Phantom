from PyQt5.QtWidgets import QFileDialog, QGroupBox, QPushButton
from PyQt5.QtGui import QFontMetrics
from PyQt5.QtCore import Qt, QThread

import pandas as pd
import os

from utils.nii_processing import save_nifti_mrs
from simulation.mrs_phantom import DigitalPhantom
from simulation.run_matlab import MatlabRunner
from simulation.basis import Basis

from gui.defaults import DEFAULT_METABOLITES
from utils.definitions import METABS, METABS_TRANSLATION
from gui.worker import Worker
from gui.widgets import SaveNiftiDialog

class Controller:
    def __init__(self, gui):
        self.gui = gui

    def select_metab_df_path(self):
        filepath, _ = QFileDialog.getOpenFileName(
            self.gui, "Select Metabolite CSV", "", "CSV Files (*.csv)"
        )
        if filepath: 

            # Load metabolites from the selected path
            self.gui.metab_df_path = filepath
            self.load_metabolites_from_path(filepath)
    
    def select_path2basis(self):
        filepath, _ = QFileDialog.getOpenFileName(
            self.gui, "Select Basis Set File", "", "Basis Files (*.BASIS)"
        )
        if filepath:
            # Set the basis set path in the GUI
            self.gui.path2basis = filepath
            self.set_elided_text(self.gui.basisset_widget.custom_basis_path_label, filepath)
            self.gui.log_message(f"✅ Basis set file selected: {filepath}")
            self.load_metabolites_from_basis(filepath)

    def set_elided_text(self, label, text, max_width=300):
        """
        Set the elided text for a label to fit within a specified width.
        """
        metrics = QFontMetrics(label.font())
        elided_text = metrics.elidedText(text, Qt.ElideMiddle, max_width)
        label.setText(elided_text)
        label.setToolTip(text)
    
    def load_metabolites_from_path(self, path):
        try:
            self.set_elided_text(self.gui.metab_file_label, f"Selected: {path}")

            # Load the CSV file and extract metabolites
            df = pd.read_csv(path)
            metabolites = df["Metabolite"].dropna().unique().tolist()

            # Only keep metabolites that are present in the default list AND in the CSV
            # But first translate the default metabolites according to the definitions
            translated_metabolites = [METABS_TRANSLATION.get(metab, metab) for metab in DEFAULT_METABOLITES]
            metabolites = [metab for metab in metabolites if metab in translated_metabolites]

            metabolites.sort(key=lambda x: x.lower())
            self.gui.update_metabolites(metabolites)
        except Exception as e:
            self.gui.log_message(f"❌ Error loading metabolites from {path}: {e}")
    
    def load_metabolites_from_basis(self, path):
        try:
            basis_settings = self.gui.basisset_widget.get_all_settings()
            metab_df = pd.read_csv(self.gui.metab_df_path)
            basis = Basis(path2basis=path, metab_df=metab_df, vendor=basis_settings['vendor'],
                          localization=basis_settings['localization'], TE=basis_settings['TE'])
            metabolites_basis = basis.names
            metab_df = pd.read_csv(self.gui.metab_df_path)
            metabolites_df = metab_df["Metabolite"].dropna().unique().tolist()

            # Only keep metabolites that are present in the default list AND in the basis set AND in the CSV
            translated_metabolites = [METABS_TRANSLATION.get(metab, metab) for metab in DEFAULT_METABOLITES]
            metabolites = [metab for metab in metabolites_basis if metab in translated_metabolites and metab in metabolites_df]
            metabolites.sort(key=lambda x: x.lower())
            self.gui.update_metabolites(metabolites)

        except Exception as e:
            self.gui.log_message(f"❌ Error loading metabolites from basis set {path}: {e}")

    def generate_basis_set(self):
        settings = self.gui.basisset_widget.get_all_settings()

        # When custom basis set is selected, use the custom path
        if self.gui.basisset_widget.source_selector.currentText() == "Load from file":
            if not self.gui.path2basis:
                self.gui.log_message("⚠️ No basis set file selected!", color="red")
                return
            return
    
        # Otherwise, generate the basis set based on the selected settings
        expected_file_name = f"LCModel_{settings['vendor']}_UnEdited_{settings['localization']}_TE{settings['TE']}.BASIS"
        
        # If the file exists, no need to generate
        if os.path.exists(os.path.join(self.gui.basis_set_dir, expected_file_name)):
            self.gui.log_message(f"📂 Basis set file already exists: {expected_file_name}")
            self.gui.log_message("No need to generate a new one.")
        else:
            # Generate the basis set
            self.gui.log_message(f"⚙️ Generating basis set: {expected_file_name}")
            self.gui.log_message("This may take a while...")
            # Call the MatlabRunner to generate the basis set
            matlab_runner = MatlabRunner()
            matlab_runner.generate_basis_set(settings)
            self.gui.log_message(f"✅ Basis set generated: {expected_file_name}")

        # Update the path2basis in the GUI
        self.gui.path2basis = os.path.join(self.gui.basis_set_dir, expected_file_name)
    
    def load_phantom_data(self):
        self.gui.log_message("🔄 Loading phantom data...")
        skeleton_name = self.gui.skeleton_dropdown.currentText()
        
        self.generate_basis_set()

        simulation_params = self.gui.simulation_settings.get_settings()
        basis_params = self.gui.basisset_widget.get_all_settings()
        simulation_params = {**simulation_params, **basis_params}

        phantom = DigitalPhantom(
            skeleton=skeleton_name,
            path2metabs=self.gui.metab_df_path,
            path2basis=self.gui.path2basis,
            simulation_params=simulation_params,
            gui=self.gui,
        )
        
        volume = phantom.get_image_data()
        labels = phantom.get_label_data().numpy()
        lipid_mask = phantom.lipid_mask.numpy()
        B0_map = phantom.get_B0_data()
        voxel_spacing = phantom.spacing

        volumes = {
            "Image": volume,
            "Labels": labels,
            "Lipid mask": lipid_mask
        }

        # If B0 map is available, add it to the volumes
        if B0_map is not None:
            volumes["B0 map"] = B0_map

        return (phantom, volumes, voxel_spacing)

    def on_phantom_loaded(self, result):
        phantom, volumes, voxel_spacing = result

        self.phantom = phantom

        self.data_volume = volumes

        self.gui.phantom_viewer.set_volume(self.data_volume, voxel_spacing)

        self.gui.log_message(f"✅ Phantom loaded with skeleton: {self.gui.skeleton_dropdown.currentText()}")
        self.gui.log_message(f"✅ Metabolite DataFrame path: {self.gui.metab_df_path}")
        self.gui.log_message(f"✅ Basis set path: {self.gui.path2basis}")

        self.disable_left_panel()
    
    def load_phantom_debug(self):
        # This method is for debugging purposes
        result = self.load_phantom_data()
        self.on_phantom_loaded(result)

    def generate_spectrum(self):
        # Get the selected volume coordinates from the GUI
        coords = self.gui.phantom_viewer.get_selected_coords()
        if coords is None:
            print("No coordinates selected.")
            return
        
        # Set simulation parameters and update the phantom
        simulation_params = self.gui.simulation_settings.get_settings()
        basis_params = self.gui.basisset_widget.get_all_settings()
        simulation_params = {**simulation_params, **basis_params}

        self.phantom.update_phantom(simulation_params)

        volume = self.phantom.get_image_data()
        labels = self.phantom.get_label_data().numpy()
        lipid_mask = self.phantom.lipid_mask.numpy()

        volumes = {
            "Image": volume,
            "Labels": labels,
            "Lipid mask": lipid_mask
        }

        # If B0 map is available, add it to the volumes
        B0_map = self.phantom.get_B0_data()
        if B0_map is not None:
            volumes["B0 map"] = B0_map

        self.data_volume = volumes

        # Update the data in the phantom viewer
        self.gui.phantom_viewer.set_volume(self.data_volume, self.phantom.spacing, refresh=False)

        # Simulate the MRS data
        self.gui.log_message("🧠 Simulating MRS data...")
        simulated_spec, components, time_axis, ppm_axis = self.phantom.simulate_data(coords)


        return (simulated_spec, components, time_axis, ppm_axis)

    def on_spectrum_generated(self, result):
        # Unpack the result
        simulated_spec, components, time_axis, ppm_axis = result

        # Plot the simulated spectrum in the GUI
        self.gui.spectrum_plot_widget.set_data(simulated_spec, ppm_axis, time_axis, components)

    def generate_spectrum_debug(self):
        # This method is for debugging purposes
        result = self.generate_spectrum()
        self.on_spectrum_generated(result)

    def reset_ui(self):
        # Enable the left panel elements
        self.gui.skeleton_panel.setEnabled(True)
        self.gui.data_selection_panel.setEnabled(True)
        self.gui.basisset_panel.setEnabled(True)

        # Reset the plots
        self.gui.phantom_viewer.reset()
        self.gui.spectrum_plot_widget.reset()

        # Clear message log
        self.gui.message_box.clear()
        self.gui.log_message("✅ Reset completed!")

    def disable_left_panel(self):
        # Disable the left panel elements
        self.gui.skeleton_panel.setEnabled(False) 
        self.gui.data_selection_panel.setEnabled(False)
        self.gui.basisset_panel.setEnabled(False)
        
        print("Left panel disabled.")
    
    def start_long_task(self, func, on_finished=None):
        self.worker = Worker(func)
        self.worker.message.connect(self.gui.log_message)
        self.worker.progress.connect(self.gui.update_progress_bar)

        # Start busy animation if no progress reported
        self.gui.update_progress_bar()  # Busy start

        if on_finished is not None:
            # Make a wrapper to first call your callback, then stop progress bar
            def finished_wrapper(result):
                on_finished(result)
                self.gui.update_progress_bar(100)  # Stop progress bar
            self.worker.finished.connect(finished_wrapper)
        else:
            # If no callback, just stop the progress bar
            self.worker.finished.connect(lambda _: self.gui.update_progress_bar(100))

        self.worker.start()
    
    def open_save_popup(self):
        spectrum_widget = self.gui.spectrum_plot_widget
        if spectrum_widget is None:
            self.gui.log_message("❌ No spectrum data available!", color="red")
            return
        
        components = spectrum_widget.components
        self.save_dialog = SaveNiftiDialog(components, self.save_selected_nifti)
        self.save_dialog.exec_()

    def save_selected_nifti(self, total_selected, individual_selected):
        components = self.gui.spectrum_plot_widget.components


        total_selection = {k: components[k] for k in total_selected}
        individual_selection = {k: components[k] for k in individual_selected}

        save_path = QFileDialog.getExistingDirectory(self.gui, "Select Directory to Save Files")
        if not save_path:
            return

        affine = self.phantom.affine
        dwelltime = self.phantom.basis_set.dwelltime
        spec_freq = self.phantom.basis_set.cf

        simulation_params = {}

        # Add vendor, localization, TE, skeleton to simulation_params
        simulation_params['skeleton'] = self.gui.skeleton_dropdown.currentText()
        simulation_params['path2metabs'] = self.gui.metab_df_path
        simulation_params['basis_set_dir'] = self.gui.basis_set_dir
        simulation_params['metabs'] = self.gui.basisset_widget.metabolite_checkboxes.get_selected_metabolites()

        simulation_params["basis_set_settings"] = {
            "vendor": self.gui.basisset_widget.vendor_dropdown.currentText(),
            "localization": self.gui.basisset_widget.localization_dropdown.currentText(),
            "TE": int(self.gui.basisset_widget.te_input.text()),
        }

        # Get simulation parameters from the GUI
        simulation_params["simulation_params"] = self.gui.simulation_settings.get_settings()

        # Extract voxel definition from the phantom
        coords = self.gui.phantom_viewer.get_selected_coords()
        voxel_size = self.gui.phantom_viewer.get_voxel_size()
        voxel_size_mm = self.gui.phantom_viewer.get_voxel_size_mm()
        voxel_definition = {
            "coords": coords,
            "size": voxel_size,
            "size_mm": voxel_size_mm
        }
        simulation_params['voxel_definitions'] = voxel_definition

        save_nifti_mrs(save_path, total_selection, individual_selection, components, affine, dwelltime, spec_freq, nucleus='1H',
                       sim_params=simulation_params)


