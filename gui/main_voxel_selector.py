from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QPushButton, QTableWidget, QTableWidgetItem, QHeaderView, QApplication
from gui.widgets import PhantomViewer, VoxelSizeWidget
import os
import json
import random
import string
import time
import copy

class PhantomWindow(QMainWindow):
    def __init__(self, phantom, config, output_dir=None):
        super().__init__()
        self.setWindowTitle("Phantom Viewer")
        
        # Store the selected voxel coordinates
        self.selected_voxels = []   # List to store selected voxels
        self.voxel_sizes = []       # List to store voxel sizes
        self.voxel_sizes_mm = []    # List to store voxel sizes in mm

        self.original_config = config  # Save a copy of the loaded config
        self.output_dir = output_dir if output_dir else os.path.dirname(__file__)

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir, exist_ok=True)

        # Create a widget to hold the PhantomViewer
        widget = QWidget(self)
        layout = QVBoxLayout()

        # Instantiate the VoxelSizeWidget and PhantomViewer
        self.voxel_size_widget = VoxelSizeWidget()
        layout.addWidget(self.voxel_size_widget)

        volume = phantom.get_image_data()
        labels = phantom.get_label_data().numpy()
        lipid_mask = phantom.lipid_mask.numpy()

        volumes = {
            "Image": volume,
            "Labels": labels,
            "Lipid mask": lipid_mask
        }

        self.viewer = PhantomViewer(volumes, voxel_spacing=phantom.spacing, voxel_size_widget=self.voxel_size_widget)
        self.viewer.set_volume(volumes, phantom.spacing)
        layout.addWidget(self.viewer)

        # Button to extract voxel coordinates
        self.extract_voxel_button = QPushButton("Add Selected Voxel")
        self.extract_voxel_button.clicked.connect(self.add_selected_voxel)
        layout.addWidget(self.extract_voxel_button)

        # Create a table to display selected voxels
        self.voxel_table = QTableWidget()
        self.voxel_table.setColumnCount(3)
        self.voxel_table.setHorizontalHeaderLabels(["Coordinates", "Size (voxels)", "Size (mm)"])
        layout.addWidget(self.voxel_table)

        self.voxel_table.horizontalHeader().setStretchLastSection(True)
        self.voxel_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.voxel_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.voxel_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)

        # Button to remove selected voxel from the table
        self.remove_voxel_button = QPushButton("Remove Selected Voxel")
        self.remove_voxel_button.clicked.connect(self.remove_selected_voxel)
        layout.addWidget(self.remove_voxel_button)

        # Button to confirm and save
        self.save_config_button = QPushButton("Confirm and Save Config")
        self.save_config_button.clicked.connect(self.save_selected_voxels)
        layout.addWidget(self.save_config_button)

        # Set the layout for the window
        widget.setLayout(layout)
        self.setCentralWidget(widget)

    def add_selected_voxel(self):
        """Adds the currently selected voxel coordinates and voxel size."""
        coords = self.viewer.get_selected_coords()
        voxel_size = self.viewer.get_voxel_size()
        voxel_size_mm = self.viewer.get_voxel_size_mm()

        if coords is not None:
            self.selected_voxels.append(coords)
            self.voxel_sizes.append(voxel_size)
            self.voxel_sizes_mm.append(voxel_size_mm)

            print(f"Added voxel: {coords} with size {voxel_size}, size in mm: {voxel_size_mm}")

            # Add a new row to the table
            row_position = self.voxel_table.rowCount()
            self.voxel_table.insertRow(row_position)

            # Fill the cells
            self.voxel_table.setItem(row_position, 0, QTableWidgetItem(str(coords)))
            self.voxel_table.setItem(row_position, 1, QTableWidgetItem(str(voxel_size)))
            self.voxel_table.setItem(row_position, 2, QTableWidgetItem(str(voxel_size_mm)))
        else:
            print("No voxel selected.")

    def remove_selected_voxel(self):
        """Removes the selected voxel from the table and from the stored lists."""
        selected_rows = self.voxel_table.selectionModel().selectedRows()
        if not selected_rows:
            print("No voxel selected for removal.")
            return

        for selected_row in sorted(selected_rows, reverse=True):  # Reverse so deletion doesn't shift indices
            row_idx = selected_row.row()

            # Remove from internal lists
            del self.selected_voxels[row_idx]
            del self.voxel_sizes[row_idx]
            del self.voxel_sizes_mm[row_idx]

            # Remove from table
            self.voxel_table.removeRow(row_idx)

        print(f"Removed {len(selected_rows)} voxel(s).")

    def save_selected_voxels(self):
        """Save the selected voxels and voxel sizes to a new config file."""
        if not self.selected_voxels:
            print("No voxels selected; nothing to save.")
            return

        print("Saving selected voxel definitions to new config file...")

        new_config = copy.deepcopy(self.original_config)

        # Define a new "voxel_definitions" field
        voxel_definitions = []
        for coord, size, size_mm in zip(self.selected_voxels, self.voxel_sizes, self.voxel_sizes_mm):
            voxel_definitions.append({
                "coords": coord,
                "size": size,
                "size_mm": size_mm
            })

        new_config['voxel_definitions'] = voxel_definitions

        # Create a unique ID (YearMonthDay-random) for folder name
        date_str = time.strftime("%Y%m%d")
        unique_id = ''.join(random.choices(string.ascii_lowercase + string.digits, k=4))
        # Create a unique folder name
        folder_name = f"simulation_{date_str}_{unique_id}"
        # Create the output directory if it doesn't exist
        os.makedirs(os.path.join(self.output_dir, folder_name), exist_ok=True)
        # Save the new config file in the output directory
        
        # Build new config file name
        new_config_path = os.path.join(self.output_dir, folder_name, f"config.json")

        with open(new_config_path, 'w') as f:
            json.dump(new_config, f, indent=4)

        print(f"Saved updated config file to: {new_config_path}")
        self.config_path = new_config_path
        self.close()
        QApplication.quit()