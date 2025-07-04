from PyQt5.QtWidgets import QVBoxLayout, QCheckBox, QPushButton, QLabel, QHBoxLayout, QGroupBox, QSizePolicy, \
                            QGridLayout, QComboBox, QSlider, QWidget, QLineEdit, QDialog, QFrame, QSpacerItem, \
                            QTabWidget, QMenu, QToolButton, QAction, QFileDialog

from PyQt5.QtGui import QDoubleValidator, QIntValidator, QFontMetrics
from PyQt5.QtCore import Qt, QRectF, QLocale
from PyQt5.QtCore import pyqtSignal

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from simulation.basis import Basis

import pyqtgraph as pg
import numpy as np
import torch

from gui.defaults import *

class MetaboliteCheckboxes:
    def __init__(self, metabolite_list, columns=4):
        self.columns = columns
        self.layout = QVBoxLayout()
        
        # Add the Select/Deselect All checkbox
        self.select_all_checkbox = QCheckBox("Select/Deselect All")
        self.select_all_checkbox.stateChanged.connect(self.toggle_all_metabolites)
        self.layout.addWidget(self.select_all_checkbox)
        
        # Create the grid layout for the metabolites
        self.grid_layout = QGridLayout()
        self.checkboxes = {}
        
        for idx, metab in enumerate(metabolite_list):
            row = idx // columns
            col = idx % columns
            checkbox = QCheckBox(metab)
            checkbox.setChecked(True)  # Default to checked
            # if metab in ("Cr", "H2O"):
            #     checkbox.setEnabled(False)  # Disable Cr and H2O
            self.checkboxes[metab] = checkbox
            self.grid_layout.addWidget(checkbox, row, col)
        
        # Add the grid layout to the main layout
        self.layout.addLayout(self.grid_layout)
    
    def toggle_all_metabolites(self, state):
        """
        This method will select/deselect all metabolites based on the state of the checkbox.
        If the checkbox is checked, it will select all metabolite checkboxes, otherwise deselect all.
        """
        for metab, checkbox in self.checkboxes.items():
            # if metab in ("Cr", "H2O"):
            #     continue  # Don't change Cr and H2O
            checkbox.setChecked(state == 2)

    def get_selected_metabolites(self):
        """
        Returns a list of selected metabolites based on the checked state of the checkboxes.
        """
        return [metab for metab, cb in self.checkboxes.items() if cb.isChecked()]
    
    def get_deselected_metabolites(self):
        """
        Returns a list of deselected metabolites based on the unchecked state of the checkboxes.
        """
        return [metab for metab, cb in self.checkboxes.items() if not cb.isChecked()]
    
    def update_metabolites(self, new_metabolites):
        """
        Updates the list of metabolites and refreshes the checkboxes.
        """
        # Clear existing checkboxes
        for checkbox in self.checkboxes.values():
            self.grid_layout.removeWidget(checkbox)
            checkbox.deleteLater()
        
        # Clear the dictionary
        self.checkboxes.clear()
        
        # Add new checkboxes
        for idx, metab in enumerate(new_metabolites):
            row = idx // self.columns
            col = idx % self.columns
            checkbox = QCheckBox(metab)
            checkbox.setChecked(True)  # Default to checked
            # if metab in ("Cr", "H2O"):
            #     checkbox.setEnabled(False)  # Disable Cr and H2O
            self.checkboxes[metab] = checkbox
            self.grid_layout.addWidget(checkbox, row, col)

class BasisSetSettingsWidget(QWidget):
    def __init__(self, default_metabolites, columns=4, parent=None):
        super().__init__(parent)
        self.layout = QVBoxLayout(self)

        # 0. Basis Set Source
        self.layout.addWidget(QLabel("Basis Set Source"))
        self.source_selector = QComboBox()
        self.source_selector.addItems(["Generate Basis Set", "Load from file"])
        self.source_selector.currentIndexChanged.connect(self.toggle_basis_input_mode)
        self.layout.addWidget(self.source_selector)

        # 1. Custom basis set load (initially hidden)
        self.custom_basis_path_label = QLabel("No custom basis set loaded")
        self.custom_basis_path_label.setVisible(False)
        self.layout.addWidget(self.custom_basis_path_label)

        self.browse_basis_button = QPushButton("Browse Basis Set File")
        self.browse_basis_button.setVisible(False)
        # self.browse_basis_button.clicked.connect(self.controller.set_path2basis)
        self.layout.addWidget(self.browse_basis_button)

        # 2. Metabolite selection      
        self.layout.addWidget(QLabel("Select Metabolites for Basis Set"))
        self.metabolite_checkboxes = MetaboliteCheckboxes(default_metabolites, columns=columns)
        self.layout.addLayout(self.metabolite_checkboxes.layout)

        # 3. Vendor selection
        self.vendor_label = QLabel("Select Vendor")
        self.layout.addWidget(self.vendor_label)
        self.vendor_dropdown = QComboBox()
        self.vendor_dropdown.addItems(["Universal_Philips", "Universal_Siemens"])
        self.layout.addWidget(self.vendor_dropdown)

        # 4. Localization method
        self.layout.addWidget(QLabel("Localization Method"))
        self.localization_dropdown = QComboBox()
        self.localization_dropdown.addItems(["PRESS", "sLASER"])
        self.layout.addWidget(self.localization_dropdown)

        # 5. Echo time (TE)
        self.layout.addWidget(QLabel("Echo Time (ms)"))
        self.te_input = QLineEdit()
        self.te_input.setText(str(DEFAULT_ECHO_TIME_MS))
        self.te_input.setValidator(QIntValidator(1, 10000))  # Optional: restrict to 1–10000 ms
        self.layout.addWidget(self.te_input)

    def toggle_basis_input_mode(self, index):
        is_file_mode = index == 1  # 0 = Generate, 1 = Load from file

        self.browse_basis_button.setVisible(is_file_mode)
        self.custom_basis_path_label.setVisible(is_file_mode)

        self.vendor_label.setVisible(not is_file_mode)
        self.vendor_dropdown.setVisible(not is_file_mode)
        self.vendor_dropdown.setEnabled(not is_file_mode)

    def get_all_settings(self):
        """Returns all the currently set parameters in a dictionary."""
        settings = {
            "selected_metabolites": self.metabolite_checkboxes.get_selected_metabolites(),
            "vendor": self.vendor_dropdown.currentText(),
            "localization": self.localization_dropdown.currentText(),
            "TE": int(self.te_input.text()),
        }
        return settings

class ImageSlider(QSlider):
    def __init__(self, axis=0, parent=None):
        super().__init__(parent)
        self.axis = axis
        self.setOrientation('horizontal')
        self.setRange(0, 100)  # Placeholder range, will be set later
        self.valueChanged.connect(self.on_value_changed)

    def on_value_changed(self, value):
        if hasattr(self.parent(), 'set_slice'):
            self.parent().set_slice(value)

class VoxelSizeWidget(QWidget):
    # Signal emitted when user finishes editing voxel size
    voxelSizeChanged = pyqtSignal(float, float, float)

    def __init__(self, parent=None):
        super().__init__(parent)

        self.x_input = QLineEdit()
        self.y_input = QLineEdit()
        self.z_input = QLineEdit()

        validator = QDoubleValidator(0.0, 10000.0, 2, parent=self)
        validator.setLocale(QLocale(QLocale.C))            

        # Set default values
        for inp in (self.x_input, self.y_input, self.z_input):
            inp.setValidator(validator)
            inp.setFixedWidth(60)

        self.x_input.setText(str(DEFAULT_VOXEL_SIZE_MM))
        self.y_input.setText(str(DEFAULT_VOXEL_SIZE_MM))
        self.z_input.setText(str(DEFAULT_VOXEL_SIZE_MM))

        # Connect editingFinished to emit signal
        for inp in (self.x_input, self.y_input, self.z_input):
            inp.editingFinished.connect(self._on_edit_finished)

        def pair(label_text, field):
            layout = QHBoxLayout()
            layout.setSpacing(5)
            layout.setContentsMargins(0,0,0,0)
            layout.addWidget(QLabel(label_text))
            layout.addWidget(field)
            return layout

        main_layout = QHBoxLayout()
        for lbl, fld in zip(["X (mm):", "Y (mm):", "Z (mm):"],
                            [self.x_input, self.y_input, self.z_input]):
            main_layout.addLayout(pair(lbl, fld))
        self.setLayout(main_layout)

    def _on_edit_finished(self):
        """
        1) Read whatever the user typed.
        2) Round each to 2 decimals.
        3) Re-write the fields (forcing two decimals, adding trailing zeroes).
        4) Emit the signal with the rounded values.
        """
        try:
            x = float(self.x_input.text())
            y = float(self.y_input.text())
            z = float(self.z_input.text())
        except ValueError:
            return  # invalid text, ignore

        # Round to 2 decimals
        x2, y2, z2 = round(x, 2), round(y, 2), round(z, 2)

        #  Block signals so setText() doesn’t re-enter us:
        for fld in (self.x_input, self.y_input, self.z_input):
            fld.blockSignals(True)

        # Format with two decimals (this will turn “10.4” → “10.40”, or “5” → “5.00”,
        # and “3.14159” → “3.14”)
        self.x_input.setText(f"{x2:.2f}")
        self.y_input.setText(f"{y2:.2f}")
        self.z_input.setText(f"{z2:.2f}")

        for fld in (self.x_input, self.y_input, self.z_input):
            fld.blockSignals(False)

        # Finally emit the normalized values
        self.voxelSizeChanged.emit(x2, y2, z2)

    def set_voxel_size(self, size_tuple):
        # Block signals to avoid recursive emit
        for inp in (self.x_input, self.y_input, self.z_input):
            inp.blockSignals(True)
        x, y, z = size_tuple
        self.x_input.setText(f"{x:.2f}")
        self.y_input.setText(f"{y:.2f}")
        self.z_input.setText(f"{z:.2f}")
        for inp in (self.x_input, self.y_input, self.z_input):
            inp.blockSignals(False)


class SimulationSettingsWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)

        int_validator = QIntValidator(1, 100000, self)
        float_validator = QDoubleValidator(0.0, 1e6, 3, self)

        def make_line_edit(validator, default_value, width=80):
            field = QLineEdit()
            field.setValidator(validator)
            field.setFixedWidth(width)
            field.setText(str(default_value))
            return field
        
        def file_picker_row(label_text, default_path="", file_filter="JSON Files (*.json)"):
            layout = QHBoxLayout()
            layout.setSpacing(10)
            layout.setContentsMargins(0, 0, 0, 0)

            label = QLabel(label_text)
            label.setFixedWidth(180)
            label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)

            path_field = QLineEdit()
            path_field.setText(default_path)
            path_field.setFixedWidth(200)

            button = QPushButton("Browse")
            button.setFixedWidth(80)

            def browse():
                file_path, _ = QFileDialog.getOpenFileName(None, "Select MM JSON File", "", file_filter)
                if file_path:
                    path_field.setText(file_path)

            button.clicked.connect(browse)

            layout.addWidget(label)
            layout.addWidget(path_field)
            layout.addWidget(button)
            layout.addStretch()

            return layout, path_field

        def labeled_row(label_text, field):
            layout = QHBoxLayout()
            layout.setSpacing(10)
            layout.setContentsMargins(0, 0, 0, 0)
            label = QLabel(label_text)
            label.setFixedWidth(180)
            label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
            layout.addWidget(label)
            layout.addWidget(field)
            layout.addStretch()
            return layout

        # New helper for min/max pairs
        def labeled_min_max_row(label_text, min_field, max_field):
            layout = QHBoxLayout()
            layout.setSpacing(10)
            layout.setContentsMargins(0, 0, 0, 0)
            label = QLabel(label_text)
            label.setFixedWidth(130)
            label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
            layout.addWidget(label)
            layout.addWidget(QLabel("Min:"))
            layout.addWidget(min_field)
            layout.addWidget(QLabel("Max:"))
            layout.addWidget(max_field)
            layout.addStretch()
            return layout

        # === Inputs for main simulation ===
        self.points_input = make_line_edit(int_validator, DEFAULT_SPECTRAL_POINTS)
        self.tr_input = make_line_edit(float_validator, DEFAULT_REPETITION_TIME_MS)
        self.bw_input = make_line_edit(float_validator, DEFAULT_BW)
        self.noise_input = make_line_edit(float_validator, DEFAULT_NOISE_LEVEL)

        # === MM inputs ===
        self.mm_level_input = make_line_edit(float_validator, DEFAULT_MM_LEVEL)
        mm_json_row, self.mm_json_file_input = file_picker_row("MM JSON File:", default_path=DEFAULT_MM_JSON_PATH)


        # === Lipid inputs ===
        self.lipid_amp_factor_input = make_line_edit(float_validator, DEFAULT_LIPID_AMP_FACTOR)
        self.lipid_sigma_input = make_line_edit(float_validator, DEFAULT_LIPID_SIGMA)
        self.lipid_phase_min_input = make_line_edit(float_validator, DEFAULT_LIPID_PHASE_MIN)
        self.lipid_phase_max_input = make_line_edit(float_validator, DEFAULT_LIPID_PHASE_MAX)
        self.lipid_lw_min_input = make_line_edit(float_validator, DEFAULT_LIPID_LW_MIN)
        self.lipid_lw_max_input = make_line_edit(float_validator, DEFAULT_LIPID_LW_MAX)

        # === Water inputs ===
        self.water_amp_factor_input = make_line_edit(float_validator, DEFAULT_WATER_AMP_FACTOR)
        self.water_phase_min_input = make_line_edit(float_validator, DEFAULT_WATER_PHASE_MIN)
        self.water_phase_max_input = make_line_edit(float_validator, DEFAULT_WATER_PHASE_MAX)
        self.water_damping_min_input = make_line_edit(float_validator, DEFAULT_WATER_DAMPING_MIN)
        self.water_damping_max_input = make_line_edit(float_validator, DEFAULT_WATER_DAMPING_MAX)

        # === Shim inputs ===
        self.shim_amplitude_hz_input = make_line_edit(float_validator, DEFAULT_SHIM_AMPLITUDE_HZ)
        self.shim_corr_length_input = make_line_edit(float_validator, DEFAULT_SHIM_CORR_LENGTH)
        self.shim_boundary_amp_factor_input = make_line_edit(float_validator, DEFAULT_SHIM_BOUNDARY_AMP_FACTOR)
        self.shim_boundary_smoothing_input = make_line_edit(float_validator, DEFAULT_SHIM_BOUNDARY_SMOOTHING)

        # === Tabs ===

        # Main tab
        simulation_tab = QWidget()
        simulation_layout = QVBoxLayout()
        simulation_layout.addLayout(labeled_row("Spectral Points:", self.points_input))
        simulation_layout.addLayout(labeled_row("TR (ms):", self.tr_input))
        simulation_layout.addLayout(labeled_row("Bandwidth (Hz):", self.bw_input))
        simulation_layout.addLayout(labeled_row("Noise Level:", self.noise_input))
        simulation_layout.addStretch()
        simulation_tab.setLayout(simulation_layout)

        # MM tab
        mm_tab = QWidget()
        mm_layout = QVBoxLayout()
        mm_layout.addLayout(labeled_row("MM Level:", self.mm_level_input))
        mm_layout.addLayout(mm_json_row)
        mm_layout.addStretch()
        mm_tab.setLayout(mm_layout)

        # Lipid tab
        lipid_tab = QWidget()
        lipid_layout = QVBoxLayout()
        lipid_layout.addLayout(labeled_row("Lipid amp. factor:", self.lipid_amp_factor_input))
        lipid_layout.addLayout(labeled_row("Lipid sigma (mm):", self.lipid_sigma_input))
        lipid_layout.addLayout(labeled_min_max_row("Lipid phase (deg):", self.lipid_phase_min_input, self.lipid_phase_max_input))
        lipid_layout.addLayout(labeled_min_max_row("Lipid LW (Hz):", self.lipid_lw_min_input, self.lipid_lw_max_input))
        lipid_layout.addStretch()
        lipid_tab.setLayout(lipid_layout)

        # Water tab
        water_tab = QWidget()
        water_layout = QVBoxLayout()
        water_layout.addLayout(labeled_row("Water amp. factor:", self.water_amp_factor_input))
        water_layout.addLayout(labeled_min_max_row("Water phase (deg):", self.water_phase_min_input, self.water_phase_max_input))
        water_layout.addLayout(labeled_min_max_row("Water damping (Hz):", self.water_damping_min_input, self.water_damping_max_input))
        water_layout.addStretch()
        water_tab.setLayout(water_layout)

        # Shim tab
        shim_tab = QWidget()
        shim_layout = QVBoxLayout()
        shim_layout.addLayout(labeled_row("Shim sigma (Hz):", self.shim_amplitude_hz_input))
        shim_layout.addLayout(labeled_row("Corr. length (mm):", self.shim_corr_length_input))
        shim_layout.addLayout(labeled_row("Boundary factor:", self.shim_boundary_amp_factor_input))
        shim_layout.addLayout(labeled_row("Boundary smoothing (mm):", self.shim_boundary_smoothing_input))
        shim_layout.addStretch()
        shim_tab.setLayout(shim_layout)

        tabs = QTabWidget()
        tabs.addTab(simulation_tab, "Main")
        tabs.addTab(mm_tab, "MM")
        tabs.addTab(lipid_tab, "Lipid")
        tabs.addTab(water_tab, "Water")
        tabs.addTab(shim_tab, "Shim Imperfection")

        main_layout = QVBoxLayout()
        main_layout.addWidget(tabs)
        self.setLayout(main_layout)

    def get_settings(self):
        return {
            "spectral_points": int(self.points_input.text()),
            "TR": float(self.tr_input.text()),
            "bandwidth": float(self.bw_input.text()),
            "noise_level": float(self.noise_input.text()),

            "mm_level": float(self.mm_level_input.text()),
            "mm_json_file": self.mm_json_file_input.text(),

            "lipid_amp_factor": float(self.lipid_amp_factor_input.text()),
            "lipid_sigma": float(self.lipid_sigma_input.text()),
            "lipid_phase_min": float(self.lipid_phase_min_input.text()),
            "lipid_phase_max": float(self.lipid_phase_max_input.text()),
            "lipid_lw_min": float(self.lipid_lw_min_input.text()),
            "lipid_lw_max": float(self.lipid_lw_max_input.text()),

            "water_amp_factor": float(self.water_amp_factor_input.text()),
            "water_phase_min": float(self.water_phase_min_input.text()),
            "water_phase_max": float(self.water_phase_max_input.text()),
            "water_damping_min": float(self.water_damping_min_input.text()),
            "water_damping_max": float(self.water_damping_max_input.text()),

            "shim_amplitude_hz": float(self.shim_amplitude_hz_input.text()),
            "shim_corr_length": float(self.shim_corr_length_input.text()),
            "shim_boundary_amp_factor": float(self.shim_boundary_amp_factor_input.text()),
            "shim_boundary_smoothing": float(self.shim_boundary_smoothing_input.text())
        }


class VOI:
    def __init__(self, size, position):
        """
        Initializes the VOI with a 3D size and position.
        :param size: Tuple of (width, height, depth)
        :param position: Tuple of (x, y, z) for the position of the ROI center.
        """
        self.size = size  # (width, height, depth)
        self.position = position  # (x, y, z) center position

        self.calc_bounds()

    def calc_size(self):
        """Updates the size of the VOI."""
        self.size = (
            self.x_max - self.x_min,
            self.y_max - self.y_min,
            self.z_max - self.z_min
        )

    def calc_position(self):
        """Calculates the position of the VOI based on its bounds."""
        x = (self.x_min + self.x_max) // 2
        y = (self.y_min + self.y_max) // 2
        z = (self.z_min + self.z_max) // 2
        
        self.position = (x, y, z)
    
    def calc_bounds(self):
        """Calculates the bounds of the VOI based on its size and position."""
        self.x_min = round(self.position[0] - self.size[0] // 2)
        self.x_max = round(self.position[0] + self.size[0] // 2)
        self.y_min = round(self.position[1] - self.size[1] // 2)
        self.y_max = round(self.position[1] + self.size[1] // 2)
        self.z_min = round(self.position[2] - self.size[2] // 2)
        self.z_max = round(self.position[2] + self.size[2] // 2)
    
    def print_bounds(self):
        """Prints the bounds of the VOI."""
        print(f"VOI Bounds: x_min={self.x_min}, x_max={self.x_max}, "
              f"y_min={self.y_min}, y_max={self.y_max}, "
              f"z_min={self.z_min}, z_max={self.z_max}")
    

class CustomImageView(pg.ImageView):
    def updateImage(self, *args, **kwargs):
        # Save current view range before the update
        view_range = self.getView().viewRange()
        super().updateImage(*args, **kwargs)
        # Restore the previous view range
        self.getView().setRange(xRange=view_range[0], yRange=view_range[1], padding=0)


class PhantomViewer(QWidget):
    def __init__(self, data_volumes, voxel_spacing=(1.0, 1.0, 1.0), voxel_size_widget=None):
        super().__init__()
        self.volumes = data_volumes
        self.current_volume = next(iter(data_volumes.values())) # Default to the first volume

        self.shape = self.current_volume.shape
        self.voxel_spacing = voxel_spacing # (x, y, z) spacing in mm

        self.voxel_size_widget = voxel_size_widget
        
        # Initialize VOI object
        voi_size_voxels = (
            round(DEFAULT_VOXEL_SIZE_MM / self.voxel_spacing[0]),
            round(DEFAULT_VOXEL_SIZE_MM // self.voxel_spacing[1]),
            round(DEFAULT_VOXEL_SIZE_MM // self.voxel_spacing[2])
        )

        self.voi = VOI(size=voi_size_voxels,
                       position=(self.shape[0] // 2, self.shape[1] // 2, self.shape[2] // 2))
        
        self._syncing_roi = False  # Flag to prevent recursive ROI updates

        self.init_ui()

    def init_ui(self):
        main_layout = QVBoxLayout(self)

        # Volume Selector - Allow the user to switch between available volumes
        self.view_selector = QComboBox()
        self.view_selector.addItems(self.volumes.keys())  # Populate with volume names
        self.view_selector.setCurrentIndex(0)
        self.view_selector.currentIndexChanged.connect(self.change_view)

        main_layout.addWidget(self.view_selector)

        views_layout = QHBoxLayout()
        self.axial_view = CustomImageView()
        self.coronal_view = CustomImageView()
        self.sagittal_view = CustomImageView()

        for view in [self.axial_view, self.coronal_view, self.sagittal_view]:
            view.ui.roiBtn.hide()
            view.ui.menuBtn.hide()

        views_layout.addWidget(self.axial_view)
        views_layout.addWidget(self.coronal_view)
        views_layout.addWidget(self.sagittal_view)

        sliders_layout = QHBoxLayout()

        # --- Axial slider + input ---
        self.axial_slider = QSlider(Qt.Horizontal)
        self.axial_slider.setMaximum(self.shape[2] - 1)
        self.axial_slider.setValue(self.shape[2] // 2)
        
        self.axial_input = QLineEdit(str(self.axial_slider.value()))
        self.axial_input.setFixedWidth(50)
        # only allow ints within the slider’s range
        self.axial_input.setValidator(QIntValidator(0, self.shape[2]-1, self))
        
        # sync slider → input
        self.axial_slider.valueChanged.connect(lambda v: self.axial_input.setText(str(v)))
        # sync input → slider
        self.axial_input.editingFinished.connect(
            lambda: self.axial_slider.setValue(int(self.axial_input.text() or 0))
        )
        # keep your existing image update on slider change
        self.axial_slider.valueChanged.connect(self.update_axial)

        sliders_layout.addWidget(QLabel("Axial:"))
        sliders_layout.addWidget(self.axial_slider)
        sliders_layout.addWidget(self.axial_input)


        # --- Coronal slider + input ---
        self.coronal_slider = QSlider(Qt.Horizontal)
        self.coronal_slider.setMaximum(self.shape[1] - 1)
        self.coronal_slider.setValue(self.shape[1] // 2)
        
        self.coronal_input = QLineEdit(str(self.coronal_slider.value()))
        self.coronal_input.setFixedWidth(50)
        self.coronal_input.setValidator(QIntValidator(0, self.shape[1]-1, self))
        
        self.coronal_slider.valueChanged.connect(lambda v: self.coronal_input.setText(str(v)))
        self.coronal_input.editingFinished.connect(
            lambda: self.coronal_slider.setValue(int(self.coronal_input.text() or 0))
        )
        self.coronal_slider.valueChanged.connect(self.update_coronal)

        sliders_layout.addWidget(QLabel("Coronal:"))
        sliders_layout.addWidget(self.coronal_slider)
        sliders_layout.addWidget(self.coronal_input)


        # --- Sagittal slider + input ---
        self.sagittal_slider = QSlider(Qt.Horizontal)
        self.sagittal_slider.setMaximum(self.shape[0] - 1)
        self.sagittal_slider.setValue(self.shape[0] // 2)
        
        self.sagittal_input = QLineEdit(str(self.sagittal_slider.value()))
        self.sagittal_input.setFixedWidth(50)
        self.sagittal_input.setValidator(QIntValidator(0, self.shape[0]-1, self))
        
        self.sagittal_slider.valueChanged.connect(lambda v: self.sagittal_input.setText(str(v)))
        self.sagittal_input.editingFinished.connect(
            lambda: self.sagittal_slider.setValue(int(self.sagittal_input.text() or 0))
        )
        self.sagittal_slider.valueChanged.connect(self.update_sagittal)

        sliders_layout.addWidget(QLabel("Sagittal:"))
        sliders_layout.addWidget(self.sagittal_slider)
        sliders_layout.addWidget(self.sagittal_input)


        # add to main layout
        main_layout.addLayout(views_layout)
        main_layout.addLayout(sliders_layout)
        self.setLayout(main_layout)

        # Add ROIs
        self.rois = []

        self.axial_roi = pg.RectROI(
                                    [self.voi.x_min, self.voi.y_min],
                                    [self.voi.size[0], self.voi.size[1]],
                                    maxBounds=QRectF(0, 0, self.shape[0], self.shape[1]),
                                    pen='r',
                                    movable=True,
                                    resizable=True,
                                    rotatable=False,
                                )
        self.rois.append(self.axial_roi)

        self.coronal_roi = pg.RectROI(
                                    [self.voi.x_min, self.voi.z_min],
                                    [self.voi.size[0], self.voi.size[2]],
                                    maxBounds=QRectF(0, 0, self.shape[0], self.shape[2]),
                                    pen='r',
                                    movable=True,
                                    resizable=True,
                                    rotatable=False,
                                )
        self.rois.append(self.coronal_roi)

        self.sagittal_roi = pg.RectROI(
                                    [self.voi.y_min, self.voi.z_min],
                                    [self.voi.size[1], self.voi.size[2]],
                                    maxBounds=QRectF(0, 0, self.shape[1], self.shape[2]),
                                    pen='r',
                                    movable=True,
                                    resizable=True,
                                    rotatable=False,
                                )
        self.rois.append(self.sagittal_roi)

        # Add ROIs to views
        self.axial_view.addItem(self.axial_roi)
        self.coronal_view.addItem(self.coronal_roi)
        self.sagittal_view.addItem(self.sagittal_roi)

        # Connect ROI change signals
        self.axial_roi.sigRegionChanged.connect(lambda: self.roi_moved('axial'))
        self.coronal_roi.sigRegionChanged.connect(lambda: self.roi_moved('coronal'))
        self.sagittal_roi.sigRegionChanged.connect(lambda: self.roi_moved('sagittal'))

        # Connect voxel size changes
        if self.voxel_size_widget:
            self.voxel_size_widget.voxelSizeChanged.connect(self.on_widget_size_changed)

        self.update_views()

    def update_views(self):
        self.update_axial()
        self.update_coronal()
        self.update_sagittal()

        # Print VOI bounds and size for debugging
        # print(f"VOI coordinates: {self.voi.x_min}, {self.voi.x_max}, {self.voi.y_min}, {self.voi.y_max}, {self.voi.z_min}, {self.voi.z_max}")
        # print(f"VOI size: {self.voi.size}") 
    
    def change_view(self):
        selected_view = self.view_selector.currentText()
        self.current_volume = self.volumes[selected_view]
        self.update_views()
        
        # Set intensity levels
        min_val, max_val = np.min(self.current_volume), np.max(self.current_volume)
        self.axial_view.setLevels(min_val, max_val)
        self.coronal_view.setLevels(min_val, max_val)
        self.sagittal_view.setLevels(min_val, max_val)

        self.update_views()

    def update_axial(self):
        z = self.axial_slider.value()
        
        # Save current view range
        view = self.axial_view.getView()
        view_range = view.viewRange() 
        # Update image
        self.axial_view.setImage(self.current_volume[:, :, z], autoLevels=False)
        # Restore view range
        view.setRange(xRange=view_range[0], yRange=view_range[1], padding=0)
        # Update ROI visibility
        self.update_roi_visibility('axial', z)

    def update_coronal(self):
        y = self.coronal_slider.value()
        
        view = self.coronal_view.getView()
        view_range = view.viewRange()
        
        self.coronal_view.setImage(self.current_volume[:, y, :], autoLevels=False)
        
        view.setRange(xRange=view_range[0], yRange=view_range[1], padding=0)
        # Update ROI visibility
        self.update_roi_visibility('coronal', y)


    def update_sagittal(self):
        x = self.sagittal_slider.value()
        
        view = self.sagittal_view.getView()
        view_range = view.viewRange()
        
        self.sagittal_view.setImage(self.current_volume[x, :, :], autoLevels=False)
        
        view.setRange(xRange=view_range[0], yRange=view_range[1], padding=0)
        # Update ROI visibility
        self.update_roi_visibility('sagittal', x)

    def update_roi_visibility(self, orientation, slice_index):
        # Show ROI only if it belongs to the current slice
        if orientation == 'axial':
            if self.voi.z_min <= slice_index <= self.voi.z_max:
                self.axial_roi.setVisible(True)
            else:
                self.axial_roi.setVisible(False)

        elif orientation == 'coronal':
            if self.voi.y_min <= slice_index <= self.voi.y_max:
                self.coronal_roi.setVisible(True)
            else:
                self.coronal_roi.setVisible(False)

        elif orientation == 'sagittal':
            if self.voi.x_min <= slice_index <= self.voi.x_max:
                self.sagittal_roi.setVisible(True)
            else:
                self.sagittal_roi.setVisible(False)

    def roi_moved(self, source):
        if self._syncing_roi:
            return  # Prevent recursion

        self._syncing_roi = True

        roi = getattr(self, f"{source}_roi")
        pos = list(roi.pos())
        size = list(roi.size())

        # Round the position and size to integers
        pos = [int(p) for p in pos]
        size = [int(s) for s in size]

        # Update VOI size and position based on the source ROI
        if source == 'axial':
            self.voi.x_min, self.voi.x_max = pos[0], pos[0] + size[0]
            self.voi.y_min, self.voi.y_max = pos[1], pos[1] + size[1]
        elif source == 'coronal':
            self.voi.x_min, self.voi.x_max = pos[0], pos[0] + size[0]
            self.voi.z_min, self.voi.z_max = pos[1], pos[1] + size[1]
        elif source == 'sagittal':
            self.voi.y_min, self.voi.y_max = pos[0], pos[0] + size[0]
            self.voi.z_min, self.voi.z_max = pos[1], pos[1] + size[1]

        # Update VOI position and size
        self.voi.calc_position()
        self.voi.calc_size()
    
        # Sync other ROIs and update views
        self.sync_rois()
        self.update_views()

        # Update the voxel size widget
        if self.voxel_size_widget:
            vx, vy, vz = self.voi.size
            x_mm, y_mm, z_mm = self.voxel_to_mm((vx, vy, vz))
            self.voxel_size_widget.set_voxel_size((x_mm, y_mm, z_mm))

        self._syncing_roi = False


    def sync_rois(self):
        # Temporarily block signals
        self.axial_roi.sigRegionChanged.disconnect()
        self.coronal_roi.sigRegionChanged.disconnect()
        self.sagittal_roi.sigRegionChanged.disconnect()

        self.axial_roi.setPos(self.voi.x_min, self.voi.y_min)
        self.axial_roi.setSize([self.voi.size[0], self.voi.size[1]])

        self.coronal_roi.setPos(self.voi.x_min, self.voi.z_min)
        self.coronal_roi.setSize([self.voi.size[0], self.voi.size[2]])

        self.sagittal_roi.setPos(self.voi.y_min, self.voi.z_min)
        self.sagittal_roi.setSize([self.voi.size[1], self.voi.size[2]])

        # Reconnect signals
        self.axial_roi.sigRegionChanged.connect(lambda: self.roi_moved('axial'))
        self.coronal_roi.sigRegionChanged.connect(lambda: self.roi_moved('coronal'))
        self.sagittal_roi.sigRegionChanged.connect(lambda: self.roi_moved('sagittal'))
    
    def set_volume(self, volumes_dict, new_voxel_spacing=None, refresh=True):
        """
        Load multiple volumes into the viewer in a fully modular way.
        
        Parameters:
            volumes_dict (dict): A dictionary with {name: 3D volume} entries.
            new_voxel_spacing (tuple): Optional new voxel spacing (dx, dy, dz).
        """
        if not volumes_dict:
            raise ValueError("The volumes_dict cannot be empty.")
    
        # Get current selected volume name
        current_name = self.view_selector.currentText()

        # Clear previous volumes and dropdown (without triggering signals)
        self.volumes = {}
        self.view_selector.blockSignals(True)
        self.view_selector.clear()
        self.view_selector.blockSignals(False)

        # Add new volumes (reversed to match orientation)
        for name, vol in volumes_dict.items():
            reversed_vol = vol[::-1, ::-1, ::-1].copy()
            self.volumes[name] = reversed_vol
            self.view_selector.addItem(name.capitalize())

        # Set the current volume to the one that was selected before
        if current_name in self.volumes:
            self.current_volume_name = current_name
            self.current_volume = self.volumes[current_name]
            self.view_selector.setCurrentText(current_name)
        else:
            # If the previous volume is not available, set to the first one
            self.current_volume_name = list(self.volumes.keys())[0]
            self.current_volume = self.volumes[self.current_volume_name]
            self.view_selector.setCurrentIndex(0)
        
        self.shape = self.current_volume.shape

        # Optional voxel spacing update
        if new_voxel_spacing is not None:
            self.voxel_spacing = new_voxel_spacing

        if refresh:
            # Update sliders
            self.axial_slider.setMaximum(self.shape[2] - 1)
            self.coronal_slider.setMaximum(self.shape[1] - 1)
            self.sagittal_slider.setMaximum(self.shape[0] - 1)

            self.axial_slider.setValue(self.shape[2] // 2)
            self.coronal_slider.setValue(self.shape[1] // 2)
            self.sagittal_slider.setValue(self.shape[0] // 2)

            # Update input validators
            self.axial_input.setValidator(QIntValidator(0, self.shape[2]-1, self))
            self.coronal_input.setValidator(QIntValidator(0, self.shape[1]-1, self))
            self.sagittal_input.setValidator(QIntValidator(0, self.shape[0]-1, self))

            # Reset VOI
            self.voi.position = (self.shape[0] // 2, self.shape[1] // 2, self.shape[2] // 2)
            self.voi.size = (
                round(DEFAULT_VOXEL_SIZE_MM / self.voxel_spacing[0]),
                round(DEFAULT_VOXEL_SIZE_MM / self.voxel_spacing[1]),
                round(DEFAULT_VOXEL_SIZE_MM / self.voxel_spacing[2])
            )
            self.voi.calc_bounds()

            # Update view ranges
            self.axial_view.getView().setRange(xRange=[0, self.shape[0]], yRange=[0, self.shape[1]], padding=0)
            self.coronal_view.getView().setRange(xRange=[0, self.shape[0]], yRange=[0, self.shape[2]], padding=0)
            self.sagittal_view.getView().setRange(xRange=[0, self.shape[1]], yRange=[0, self.shape[2]], padding=0)

            # Update ROIs
            self.axial_roi.setSize([self.voi.size[0], self.voi.size[1]])
            self.axial_roi.setPos([self.voi.x_min, self.voi.y_min])
            self.axial_roi.maxBounds = QRectF(0, 0, self.shape[0], self.shape[1])

            self.coronal_roi.setSize([self.voi.size[0], self.voi.size[2]])
            self.coronal_roi.setPos([self.voi.x_min, self.voi.z_min])
            self.coronal_roi.maxBounds = QRectF(0, 0, self.shape[0], self.shape[2])

            self.sagittal_roi.setSize([self.voi.size[1], self.voi.size[2]])
            self.sagittal_roi.setPos([self.voi.y_min, self.voi.z_min])
            self.sagittal_roi.maxBounds = QRectF(0, 0, self.shape[1], self.shape[2])

            # Center selection
            self.selected_voxel = (
                self.shape[0] // 2,
                self.shape[1] // 2,
                self.shape[2] // 2
            )
        
        # Set intensity levels
        min_val, max_val = np.min(self.current_volume), np.max(self.current_volume)
        self.axial_view.setLevels(min_val, max_val)
        self.coronal_view.setLevels(min_val, max_val)
        self.sagittal_view.setLevels(min_val, max_val)

        # Update all views
        self.update_views()


    def on_widget_size_changed(self, x_mm, y_mm, z_mm):
        # 1) Convert the user’s mm entry → integer voxels
        vx, vy, vz = self.mm_to_voxel((x_mm, y_mm, z_mm))

        # 2) Convert *those* voxels back → the “true” mm size
        #    (this will be snapped to whatever integer voxels you actually got)
        x_corr, y_corr, z_corr = self.voxel_to_mm((vx, vy, vz))

        # 3) Push the corrected mm back into the widget (it’ll format to two decimals)
        #    This updates the line-edits so the user immediately sees the snapped value.
        self.voxel_size_widget.set_voxel_size((x_corr, y_corr, z_corr))

        # 4) Now update your VOI in voxels and refresh the ROIs/views
        self.voi.x_min, self.voi.x_max = self.voi.x_min, self.voi.x_min + vx
        self.voi.y_min, self.voi.y_max = self.voi.y_min, self.voi.y_min + vy
        self.voi.z_min, self.voi.z_max = self.voi.z_min, self.voi.z_min + vz

        # Update the VOI size and position
        self.voi.calc_size()
        self.voi.calc_position()

        self.sync_rois()
        self.update_views()
        
    def mm_to_voxel(self, mm_tuple):
        """
        mm_tuple: (x_mm, y_mm, z_mm)
        self.voxel_spacing: (sx, sy, sz) in mm per voxel
        returns (vx, vy, vz) as ints
        """
        return tuple(
            int(round(mm_dim / sp))
            for mm_dim, sp in zip(mm_tuple, self.voxel_spacing)
        )

    def voxel_to_mm(self, vox_tuple):
        """
        vox_tuple: (vx, vy, vz)
        self.voxel_spacing: (sx, sy, sz)
        returns (x_mm, y_mm, z_mm) floats
        """
        return tuple(
            vox_dim * sp
            for vox_dim, sp in zip(vox_tuple, self.voxel_spacing)
    )

    def get_selected_coords(self):
        """
        Returns the selected coordinates of the VOI in the format:
        (x_min, x_max, y_min, y_max, z_min, z_max) based on the original volume. (not the flipped one)
        """
        # Inverting the flipped coordinates to get the original volume coordinates
        x_min = self.shape[0] - self.voi.x_max
        x_max = self.shape[0] - self.voi.x_min
        y_min = self.shape[1] - self.voi.y_max
        y_max = self.shape[1] - self.voi.y_min
        z_min = self.shape[2] - self.voi.z_max
        z_max = self.shape[2] - self.voi.z_min

        return (
            x_min, x_max,
            y_min, y_max,
            z_min, z_max
        )
    
    def get_voxel_size(self):
        """
        Returns the size of the VOI in voxels.
        """
        return self.voi.size
    
    def get_voxel_size_mm(self):
        """
        Returns the size of the VOI in mm.
        """
        return self.voxel_to_mm(self.voi.size)
    
    def reset(self):
        """
        Resets the complete plot window and removes the loaded data.
        """
        # Clear the views
        self.axial_view.clear()
        self.coronal_view.clear()
        self.sagittal_view.clear()

        # Clear the view selector
        self.view_selector.blockSignals(True)
        self.view_selector.clear()

        # Clear the volumes
        self.volumes = {"Image": np.zeros((128,128,128))}
        self.view_selector.addItem("Image")
        self.view_selector.setCurrentIndex(0)
        self.current_volume_name = list(self.volumes.keys())[0]
        self.current_volume = self.volumes[self.current_volume_name]
        self.shape = self.current_volume.shape

        # Reset the sliders and inputs
        self.axial_slider.setValue(self.shape[2] // 2)
        self.coronal_slider.setValue(self.shape[1] // 2)
        self.sagittal_slider.setValue(self.shape[0] // 2)
        self.axial_slider.setMaximum(self.shape[2] - 1)
        self.coronal_slider.setMaximum(self.shape[1] - 1)
        self.sagittal_slider.setMaximum(self.shape[0] - 1)

        self.axial_input.setText(str(self.axial_slider.value()))
        self.coronal_input.setText(str(self.coronal_slider.value()))
        self.sagittal_input.setText(str(self.sagittal_slider.value()))

        self.view_selector.blockSignals(False)
        

class MRSSpectrumWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("MRS Spectrum Viewer (pyqtgraph)")

        # Layouts
        main_layout = QVBoxLayout(self)

        # Top bar for view selection
        top_bar = QHBoxLayout()
        top_bar.addWidget(QLabel("View:"))
        self.view_selector = QComboBox()
        self.view_selector.addItems(["Spectrum Real", "Spectrum Imag", "FID Real", "FID Imag"])
        self.view_selector.currentIndexChanged.connect(self._refresh_plot)
        top_bar.addWidget(self.view_selector)
        top_bar.addStretch()

        # Plot widget
        self.plot_widget = pg.PlotWidget()
        self.plot_widget.setLabel('bottom', 'ppm')
        self.plot_widget.setLabel('left', 'Signal')
        self.plot_widget.showGrid(x=True, y=True)
        self.plot_widget.invertX(True)

        # Checkbox layout below the plot
        checkbox_frame = QFrame()
        self.checkbox_layout = QHBoxLayout()
        checkbox_frame.setLayout(self.checkbox_layout)

        # Assemble full layout
        main_layout.addLayout(top_bar)
        main_layout.addWidget(self.plot_widget)
        main_layout.addWidget(checkbox_frame)

        # Initialize storage
        self.lines = {}
        self.components = {}
        self.ppm = None
        self.spectrum = None

        # Button for component inclusion toggling
        self.include_components_button = QToolButton()
        self.include_components_button.setText("Components to Include")
        self.include_components_button.setPopupMode(QToolButton.InstantPopup)
        self.include_components_menu = QMenu()

        self.include_components_button.setMenu(self.include_components_menu)
        top_bar.addWidget(self.include_components_button)

        # Store included components set
        self.included_components = set()

    def set_data(self, spectrum, ppm_axis, time_axis, components):
        self.spectrum = spectrum.detach().cpu()
        self.ppm = ppm_axis
        self.time = time_axis
        self.components = {k: v.detach().cpu() for k, v in components.items()}

        # Initialize included components to all component keys
        self.included_components = set(self.components.keys())

        # Setup menu items
        self.include_components_menu.clear()
        for comp_name in self.components.keys():
            action = QAction(comp_name, self.include_components_menu)
            action.setCheckable(True)
            action.setChecked(True)
            action.toggled.connect(self._component_inclusion_changed)
            self.include_components_menu.addAction(action)

        self._refresh_plot()

    def _component_inclusion_changed(self, checked):
        action = self.sender()
        comp_name = action.text()
        if checked:
            self.included_components.add(comp_name)
        else:
            self.included_components.discard(comp_name)
        self._refresh_plot()

    def _refresh_plot(self):
        if self.spectrum is None or self.ppm is None:
            return

        view_mode = self.view_selector.currentText().lower()

        if "fid" in view_mode:
            signal = torch.fft.ifft(self.spectrum)
            x_axis = self.time
            self.plot_widget.setLabel('bottom', 'Time [a.u.]')
            self.plot_widget.invertX(False)
        else:
            signal = self.spectrum
            x_axis = self.ppm
            self.plot_widget.setLabel('bottom', 'ppm')
            self.plot_widget.invertX(True)

        # Compose total signal by including only selected components
        if self.components:
            total_signal = torch.zeros_like(signal)
            for comp_name, comp_data in self.components.items():
                if comp_name in self.included_components:
                    comp_signal = torch.fft.ifft(comp_data) if "fid" in view_mode else comp_data
                    total_signal += comp_signal
        else:
            total_signal = signal

        # If you want to apply the water_checkbox for legacy behavior:
        # You can disable/remove it now, or make it control the inclusion_components set.

        # Show total line
        self.plot_widget.clear()
        self.lines.clear()

        part = torch.real if "real" in view_mode else torch.imag
        y_total = part(total_signal).numpy()
        total_line = self.plot_widget.plot(x_axis, y_total, pen=pg.mkPen('w', width=2), name="Total")
        self.lines["Total"] = total_line

        # Plot components as before
        colors = ['r', 'g', 'y', 'b', 'm', 'c', 'orange', 'purple']
        for i, (name, data) in enumerate(self.components.items()):
            data_signal = torch.fft.ifft(data) if "fid" in view_mode else data
            y = part(data_signal).numpy()
            color = colors[i % len(colors)]
            comp_line = self.plot_widget.plot(x_axis, y, pen=pg.mkPen(color, width=2), name=name)
            comp_line.setVisible(False)
            self.lines[name] = comp_line

        if "spectrum" in view_mode:
            self.plot_widget.setXRange(5.0, 0.0, padding=0)
        else:
            self.plot_widget.enableAutoRange()

        self._update_checkboxes()

    def _update_checkboxes(self):
        # Clear existing checkboxes
        while self.checkbox_layout.count():
            child = self.checkbox_layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()

        # Add checkboxes: Total is checked, others are not
        for name in self.lines:
            checkbox = QCheckBox(name)
            checkbox.setChecked(name == "Total")
            checkbox.stateChanged.connect(self._toggle_line_visibility)
            self.checkbox_layout.addWidget(checkbox)

        self.checkbox_layout.addStretch()

    def _toggle_line_visibility(self, state):
        checkbox = self.sender()
        name = checkbox.text()
        visible = state == 2  # Qt.Checked
        if name in self.lines:
            self.lines[name].setVisible(visible)

    def reset(self):
        """
        Resets the complete plot window and removes the loaded data.
        """
        self.plot_widget.clear()
        self.lines.clear()
        self.components.clear()

        while self.checkbox_layout.count():
            child = self.checkbox_layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()\


class SaveNiftiDialog(QDialog):
    def __init__(self, components, save_callback):
        super().__init__()
        self.components = components
        self.save_callback = save_callback

        self.setWindowTitle("Save NIfTI-MRS Components")
        self.resize(400, 600)

        layout = QVBoxLayout(self)

        # --- TOTAL group ---
        total_group = QGroupBox("TOTAL signal: Components to include")
        total_layout = QVBoxLayout()
        self.total_checkboxes = {}
        for name in components.keys():
            cb = QCheckBox(name)
            cb.setChecked(True)
            cb.stateChanged.connect(self.update_summary)
            total_layout.addWidget(cb)
            self.total_checkboxes[name] = cb
        # Add select/deselect buttons
        total_buttons = QHBoxLayout()
        total_select_all = QPushButton("Select All")
        total_deselect_all = QPushButton("Deselect All")
        total_select_all.clicked.connect(lambda: self.set_all(self.total_checkboxes, True))
        total_deselect_all.clicked.connect(lambda: self.set_all(self.total_checkboxes, False))
        total_buttons.addWidget(total_select_all)
        total_buttons.addWidget(total_deselect_all)
        total_layout.addLayout(total_buttons)
        total_group.setLayout(total_layout)
        layout.addWidget(total_group)

        layout.addSpacerItem(QSpacerItem(0, 10, QSizePolicy.Minimum, QSizePolicy.Fixed))

        # --- INDIVIDUAL group ---
        components_group = QGroupBox("Components to save individually")
        components_layout = QVBoxLayout()
        self.components_checkboxes = {}
        for name in components.keys():
            cb = QCheckBox(name)
            cb.setChecked(False)
            cb.stateChanged.connect(self.update_summary)
            components_layout.addWidget(cb)
            self.components_checkboxes[name] = cb
        # Add select/deselect buttons
        components_buttons = QHBoxLayout()
        components_select_all = QPushButton("Select All")
        components_deselect_all = QPushButton("Deselect All")
        components_select_all.clicked.connect(lambda: self.set_all(self.components_checkboxes, True))
        components_deselect_all.clicked.connect(lambda: self.set_all(self.components_checkboxes, False))
        components_buttons.addWidget(components_select_all)
        components_buttons.addWidget(components_deselect_all)
        components_layout.addLayout(components_buttons)
        components_group.setLayout(components_layout)
        layout.addWidget(components_group)

        layout.addSpacerItem(QSpacerItem(0, 10, QSizePolicy.Minimum, QSizePolicy.Fixed))

        # --- Summary ---
        self.summary_label = QLabel("")
        self.summary_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(self.summary_label)

        # --- Save button ---
        self.save_button = QPushButton("Save Selected Components")
        self.save_button.clicked.connect(self.save_clicked)
        self.save_button.setEnabled(False)  # Initially disabled
        layout.addWidget(self.save_button)

        # Update summary for first time
        self.update_summary()

    def set_all(self, checkbox_dict, state: bool):
        for cb in checkbox_dict.values():
            cb.setChecked(state)

    def update_summary(self):
        total_selected = sum(cb.isChecked() for cb in self.total_checkboxes.values())
        individual_selected = sum(cb.isChecked() for cb in self.components_checkboxes.values())

        self.summary_label.setText(
            f"Selected for TOTAL: {total_selected} | Selected INDIVIDUALLY: {individual_selected}"
        )

        # Enable Save button only if at least one is selected
        if total_selected > 0 or individual_selected > 0:
            self.save_button.setEnabled(True)
        else:
            self.save_button.setEnabled(False)

    def save_clicked(self):
        total_selected = [name for name, cb in self.total_checkboxes.items() if cb.isChecked()]
        individual_selected = [name for name, cb in self.components_checkboxes.items() if cb.isChecked()]
        self.save_callback(total_selected, individual_selected)
        self.accept()