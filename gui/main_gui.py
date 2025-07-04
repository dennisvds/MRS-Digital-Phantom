from PyQt5.QtWidgets import QMainWindow, QVBoxLayout, QHBoxLayout, QWidget, QLabel, QPushButton, QGroupBox, QGridLayout, QComboBox
from PyQt5.QtWidgets import QGroupBox, QVBoxLayout, QLabel, QComboBox, QTextEdit, QApplication, QProgressBar
from PyQt5.QtWidgets import QSizePolicy

from gui.widgets import VoxelSizeWidget, SimulationSettingsWidget, BasisSetSettingsWidget, PhantomViewer, MRSSpectrumWidget
from gui.controller import Controller
from gui.defaults import *

import numpy as np

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        ### Initialize the main window ###
        self.setWindowTitle("MRS Phantom")
        self.setGeometry(100, 100, 800, 800)

        ### Initialize the controller ###
        self.controller = Controller(self)

        ### Variables ###
        # Paths
        self.metab_df_path = DEFAULT_METAB_DF_PATH
        self.skeleton_dir = DEFAULT_SKELETON_DIR
        self.basis_set_dir = DEFAULT_BASIS_SET_DIR
        self.path2basis = None

        self.metabs = DEFAULT_METABOLITES
        self.skeleton = DEFAULT_SKELETON
        self.echo_time_ms = DEFAULT_ECHO_TIME_MS

        self.phantom = None

        ### Create the main layout ###
        main_layout = QHBoxLayout()

        ### Create the left panel for inputs ###
        left_layout = QVBoxLayout()

        ### Create Skeleton Settings ###
        self.skeleton_panel = self.init_skeleton_settings()
        left_layout.addWidget(self.skeleton_panel)

        ### Create Data selection panel ###
        self.data_selection_panel = self.init_data_selection_panel()
        left_layout.addWidget(self.data_selection_panel)

        ### Create Basisset settings ###
        self.basisset_panel = self.init_basisset_settings()
        left_layout.addWidget(self.basisset_panel)

        ### Create buttons for functionality ###
        self.load_phantom_button = QPushButton("Load Phantom")
        self.load_phantom_button.clicked.connect(
            lambda: self.controller.start_long_task(self.controller.load_phantom_data, self.controller.on_phantom_loaded)
            )
        # # Debug button:
        # self.load_phantom_button.clicked.connect(self.controller.load_phantom_debug)

        left_layout.addWidget(self.load_phantom_button)

        self.reset_ui_button = QPushButton("Reset")
        self.reset_ui_button.clicked.connect(self.controller.reset_ui)
        left_layout.addWidget(self.reset_ui_button)

        ### Create the right panel for plots and settings ###
        right_layout = QVBoxLayout()

        ### Create the simulation settings panel ###
        simulation_panel = self.init_simulation_settings()

        ### Create the top right panel for plots ###
        top_panel = self.init_plots()
        right_layout.addWidget(top_panel, stretch=6)

        ### Create the bottom right panel ###
        bottom_panel = QHBoxLayout()

        ### Create the spectrum plot panel ###
        spectrum_panel = self.init_spectrum_plot()

        bottom_panel.addWidget(simulation_panel, stretch=2)
        bottom_panel.addLayout(spectrum_panel, stretch=3)

        right_layout.addLayout(bottom_panel, stretch=1)

        # Add panels to the main layout
        main_layout.addLayout(left_layout, stretch=1)
        main_layout.addLayout(right_layout, stretch=5)

        # Create a central widget and set the main layout
        central_widget = QWidget()
        central_widget.setLayout(main_layout)
        self.setCentralWidget(central_widget)


    def init_skeleton_settings(self):
        # Create the left panel for inputs
        skeleton_panel = QGroupBox("Skeleton Settings")
        skeleton_layout = QVBoxLayout()
        skeleton_panel.setLayout(skeleton_layout)

        # Add dropdown menu for skeleton selection
        self.skeleton_dropdown = QComboBox()
        self.skeleton_dropdown.addItems(["BigBrainMR", "MRiLab"])
        self.skeleton_dropdown.setCurrentText(DEFAULT_SKELETON)
        skeleton_layout.addWidget(self.skeleton_dropdown)
    
        return skeleton_panel
    
    def init_data_selection_panel(self):
        data_selection_panel = QGroupBox("Data Selection")
        data_selection_layout = QVBoxLayout()
        data_selection_panel.setLayout(data_selection_layout)

        data_selection_layout.addWidget(QLabel("Select Metabolite Dataframe"))
        self.metab_file_label = QLabel("No file selected")
        data_selection_layout.addWidget(self.metab_file_label)

        self.metab_df_button = QPushButton("Browse")
        self.metab_df_button.clicked.connect(self.controller.select_metab_df_path)
        data_selection_layout.addWidget(self.metab_df_button)

        return data_selection_panel
    
    def init_basisset_settings(self):
        # Create the BasisSetSettingsWidget
        basis_set_settings_panel = QGroupBox("Basis Set Settings")
        basis_set_settings_layout = QVBoxLayout()
        basis_set_settings_panel.setLayout(basis_set_settings_layout)

        self.basisset_widget = BasisSetSettingsWidget(default_metabolites=DEFAULT_METABOLITES, columns=4)
        self.basisset_widget.browse_basis_button.clicked.connect(self.controller.select_path2basis)
        basis_set_settings_layout.addWidget(self.basisset_widget)

        # Update metabolite checkboxes based on the selected metabolite dataframe
        self.controller.load_metabolites_from_path(self.metab_df_path)

        return basis_set_settings_panel

    def init_plots(self):
        # Create Axial, Coronal, and Sagittal plots of the phantom
        plot_panel = QGroupBox("Phantom Plots")
        plot_layout = QGridLayout()
        plot_panel.setLayout(plot_layout)
        plot_panel.setMinimumHeight(300)  # Set a minimum height for the group

        # Create the PhantomViewer widget
        dummy_data = {"Image": np.zeros((128, 128, 128))}
        self.phantom_viewer = PhantomViewer(dummy_data, voxel_size_widget=self.voxel_size_widget)
        plot_layout.addWidget(self.phantom_viewer)

        return plot_panel
    
    def init_simulation_settings(self):
        # Main container widget
        settings_panel = QWidget()
        main_layout = QVBoxLayout(settings_panel)

        # === Voxel Settings Group (on top) ===
        voxel_group = QGroupBox("Voxel Settings")
        voxel_layout = QHBoxLayout()
        voxel_group.setMaximumHeight(100)  # Set a maximum height for the group
        voxel_group.setLayout(voxel_layout)

        self.voxel_size_widget = VoxelSizeWidget()
        voxel_layout.addWidget(self.voxel_size_widget)

        main_layout.addWidget(voxel_group)

        # === Simulation Settings Group (below) ===
        simulation_group = QGroupBox("Simulation Settings")
        simulation_layout = QVBoxLayout()
        simulation_group.setLayout(simulation_layout)

        self.simulation_settings = SimulationSettingsWidget()
        simulation_layout.addWidget(self.simulation_settings)

        main_layout.addWidget(simulation_group)

        # === Generate Button ===
        self.generate_spectrum_button = QPushButton("Generate Spectrum")
        self.generate_spectrum_button.clicked.connect(
            lambda: self.controller.start_long_task(self.controller.generate_spectrum, self.controller.on_spectrum_generated)
            )
        # # Debug button:
        # self.generate_spectrum_button.clicked.connect(self.controller.generate_spectrum_debug)
        main_layout.addWidget(self.generate_spectrum_button)

        # === Save Button ===
        self.save_nifti_button = QPushButton("Save NIfTI-MRS")
        self.save_nifti_button.clicked.connect(self.controller.open_save_popup)
        main_layout.addWidget(self.save_nifti_button)

        # === Progress Bar === #
        self.progress_bar = QProgressBar()
        self.progress_bar.setValue(0)

        # Fix the height so layout stays consistent
        self.progress_bar.setMinimumHeight(15)
        self.progress_bar.setMaximumHeight(15)

        self.progress_bar.hide()  # Keep the space reserved, but visually hidden
        main_layout.addWidget(self.progress_bar)

        # === Message Box ===
        self.init_message_box()
        main_layout.addWidget(self.message_box)

        return settings_panel

    def init_spectrum_plot(self):
        # Create the spectrum plot panel
        spectrum_layout = QVBoxLayout()

        # Create the spectrum plot widget
        self.spectrum_plot_widget = MRSSpectrumWidget()
        spectrum_layout.addWidget(self.spectrum_plot_widget)

        return spectrum_layout
    
    def init_message_box(self):
        # Create the message box panel
        self.message_box = QTextEdit()
        self.message_box.setReadOnly(True)
        self.message_box.setPlaceholderText("Messages will appear here...")
    
    def log_message(self, text, color="black"):
        print(text)  # Optional: keep for terminal output during development
        colored_text = f'<span style="color:{color};">{text}</span>'
        self.message_box.append(colored_text)
        QApplication.processEvents()  # <-- Forces the GUI to update immediately
    
    def update_progress_bar(self, value=None):
        if not self.progress_bar.isVisible():
            self.progress_bar.show()
        if value is None:
            self.progress_bar.setRange(0, 0)  # Infinite/busy
        else:
            self.progress_bar.setRange(0, 100)
            self.progress_bar.setValue(value)
            if value >= 100:
                self.progress_bar.hide()

    
    def update_metabolites(self, metabolites):
        try:
            # Update the metabolite checkboxes with the new list
            self.basisset_widget.metabolite_checkboxes.update_metabolites(metabolites)
        except Exception as e:
            print(f"Error loading metabolites from path: {e}")
    

