####################################################################################################
#                                        visualizations.py                                         #
####################################################################################################
#                                                                                                  #
# Authors: D.M.J. van de Sande (d.m.j.v.d.sande@tue.nl)                                            #
#                                                                                                  #
# Created: 03/02/23                                                                                #
#                                                                                                  #
# Purpose: Defines methods to visualize data created using the Digital MRS Phantom.                #
#                                                                                                  #
####################################################################################################



#*************#
#   imports   #
#*************#
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import matplotlib.gridspec as gridspec
import numpy as np
import matplotlib.pyplot as plt

# Plot metabolite map in the axial, coronal, or sagittal plane
# And interactively select the slice indices to plot
def plot_metabolite_map(metab_map, metab_name, cmap='inferno', vmin=None, vmax=None, slice_idxs=None):
    """
    Plot a metabolite map. If slice indices are provided, plot a stationary image.
    Otherwise, create an interactive plot with sliders and voxel selection.

    Parameters:
    - metab_map (np.ndarray): The metabolite map.
    - metab_name (str): The name of the metabolite.
    - cmap (str): The colormap to use.
    - vmin (float): The minimum value of the colormap.
    - vmax (float): The maximum value of the colormap.
    - slice_idxs (list or None): If provided, plot a static image at these slice indices.
    """
    
    if slice_idxs:
        fig, axs = plt.subplots(1, 3, figsize=(15, 5))  # 1 row, 3 columns
        
        # Display the specified slices with a 90-degree rotation
        axs[0].imshow(np.rot90(metab_map[:, :, slice_idxs[2]]), cmap=cmap, vmin=vmin, vmax=vmax)
        axs[0].set_title(f"{metab_name} in axial plane", color='white', fontweight='bold')

        axs[1].imshow(np.rot90(metab_map[:, slice_idxs[1], :]), cmap=cmap, vmin=vmin, vmax=vmax)
        axs[1].set_title(f"{metab_name} in coronal plane", color='white', fontweight='bold')

        axs[2].imshow(np.rot90(metab_map[slice_idxs[0], :, :]), cmap=cmap, vmin=vmin, vmax=vmax)
        axs[2].set_title(f"{metab_name} in sagittal plane", color='white', fontweight='bold')

        # Add colorbar with white text
        cbar = fig.colorbar(im_axial, ax=axs, orientation='vertical', pad=0.02, fraction=0.05)
        cbar.set_label('[mM / IU]', fontsize=14, weight='bold', color='white')
        cbar.ax.tick_params(axis='y', colors='white', labelsize=12)  

    else:
        # Initial slice indices (middle slices if not provided)
        slice_idxs = [metab_map.shape[0] // 2, metab_map.shape[1] // 2, metab_map.shape[2] // 2]
        selected_voxel = None  # Stores the current selected voxel coordinates
        marker = [None, None, None]  # To store markers in each axis

        # Create figure and axes
        fig, axs = plt.subplots(1, 3, figsize=(15, 5))  # 1 row, 3 columns
        plt.subplots_adjust(bottom=0.2)  # Adjust space for sliders

        # Display the initial rotated images
        im_axial = axs[0].imshow(np.rot90(metab_map[:, :, slice_idxs[2]]), cmap=cmap, vmin=vmin, vmax=vmax)
        axs[0].set_title(f"{metab_name} in axial plane", color='white', fontweight='bold')

        im_coronal = axs[1].imshow(np.rot90(metab_map[:, slice_idxs[1], :]), cmap=cmap, vmin=vmin, vmax=vmax)
        axs[1].set_title(f"{metab_name} in coronal plane", color='white', fontweight='bold')

        im_sagittal = axs[2].imshow(np.rot90(metab_map[slice_idxs[0], :, :]), cmap=cmap, vmin=vmin, vmax=vmax)
        axs[2].set_title(f"{metab_name} in sagittal plane", color='white', fontweight='bold')

        # Adjust spacing to avoid overlap between colorbar and sliders
        plt.subplots_adjust(left=0.1, right=0.9, bottom=0.15, top=0.9, hspace=0.2)

        # Add colorbar
        cbar = fig.colorbar(im_axial, ax=axs, orientation='vertical', pad=0.02, fraction=0.05)
        cbar.set_label('[mM / IU]', fontsize=14, weight='bold', color='white')
        cbar.ax.tick_params(axis='y', colors='white', labelsize=12)  

        # Remove axes
        for ax in axs:
            ax.axis('off')

        # Create slider axes for each plane
        axcolor = 'black'  # Black background for sliders
        ax_slider_axial = plt.axes([0.1, 0.1, 0.2, 0.03], facecolor=axcolor)
        ax_slider_coronal = plt.axes([0.4, 0.1, 0.2, 0.03], facecolor=axcolor)
        ax_slider_sagittal = plt.axes([0.7, 0.1, 0.2, 0.03], facecolor=axcolor)

        # Create sliders with white bold labels
        slider_axial = Slider(ax_slider_axial, 'Axial', 0, metab_map.shape[2] - 1, valinit=slice_idxs[2], valstep=1)
        slider_coronal = Slider(ax_slider_coronal, 'Coronal', 0, metab_map.shape[1] - 1, valinit=slice_idxs[1], valstep=1)
        slider_sagittal = Slider(ax_slider_sagittal, 'Sagittal', 0, metab_map.shape[0] - 1, valinit=slice_idxs[0], valstep=1)

        # Set slider label and tick colors to white and bold
        for slider in [slider_axial, slider_coronal, slider_sagittal]:
            slider.label.set_color('white')
            slider.label.set_fontweight('bold')
            slider.valtext.set_color('white')
            slider.valtext.set_fontweight('bold')

        # Function to update the displayed images
        def update(val):
            slice_axial = int(slider_axial.val)
            slice_coronal = int(slider_coronal.val)
            slice_sagittal = int(slider_sagittal.val)

            # Update the metabolite map images (with rotation)
            im_axial.set_data(np.rot90(metab_map[:, :, slice_axial]))
            im_coronal.set_data(np.rot90(metab_map[:, slice_coronal, :]))
            im_sagittal.set_data(np.rot90(metab_map[slice_sagittal, :, :]))

            # Redraw the figure
            fig.canvas.draw_idle()

            # Clear previous markers
            for i in range(3):
                if marker[i] is not None:
                    marker[i].remove()
                    marker[i] = None

            # Handle marker display based on the selected voxel
            if selected_voxel:
                x, y, z = selected_voxel
                # Draw the marker on the correct plane
                if z == slice_axial:
                    marker[0], = axs[0].plot(x, metab_map.shape[1] - y, 'o', markersize=16, markeredgewidth=4, markeredgecolor='green', markerfacecolor='none')
                if y == slice_coronal:
                    marker[1], = axs[1].plot(x, metab_map.shape[2] - z, 'o', markersize=16, markeredgewidth=4, markeredgecolor='green', markerfacecolor='none')
                if x == slice_sagittal:
                    marker[2], = axs[2].plot(y, metab_map.shape[2] - z, 'o', markersize=16, markeredgewidth=4, markeredgecolor='green', markerfacecolor='none')

        # Attach the update function to each slider
        slider_axial.on_changed(update)
        slider_coronal.on_changed(update)
        slider_sagittal.on_changed(update)

        def on_click(event):
            nonlocal selected_voxel
            if event.inaxes not in axs:
                return

            # Determine which plane was clicked and calculate voxel coordinates
            if event.inaxes == axs[0]:  # Axial plane
                x = int(event.xdata)
                y = metab_map.shape[1] - int(event.ydata)
                z = int(slider_axial.val)  # Current axial slice
            elif event.inaxes == axs[1]:  # Coronal plane
                x = int(event.xdata)
                z = metab_map.shape[2] - int(event.ydata)
                y = int(slider_coronal.val)  # Current coronal slice
            elif event.inaxes == axs[2]:  # Sagittal plane
                y = int(event.xdata)
                z = metab_map.shape[2] - int(event.ydata)
                x = int(slider_sagittal.val)  # Current sagittal slice
            else:
                return

            # Update the selected voxel
            selected_voxel = (x, y, z)

            # Set sliders to the new selected voxel's slices
            slider_axial.set_val(z)       # Update axial slice
            slider_coronal.set_val(y)     # Update coronal slice
            slider_sagittal.set_val(x)    # Update sagittal slice



        # Connect the click event to the figure
        fig.canvas.mpl_connect('button_press_event', on_click)

        # Black background with bold white text
        fig.set_facecolor('black')
        for ax in axs:
            ax.set_facecolor('black')

    plt.show()

    return fig, axs

def plot_mrsi(spectra, metab_map, ppm_axis, metab_name='NAA', slice_idx=None, start_x=None, start_y=None):
    """
    Plot an MRSI dataset interactively. The user can navigate slices with a slider or arrow keys
    and click on the metabolite map to view spectra at specific locations.
    """
    # Create the figure and the two subplots
    fig = plt.figure(figsize=(10, 6))
    gs = gridspec.GridSpec(2, 2, height_ratios=[9, 1], width_ratios=[1.1, 1], wspace=0.6, hspace=0.4)
    ax_map = fig.add_subplot(gs[0, 0])
    ax_spec = fig.add_subplot(gs[0, 1])
    ax_slider = fig.add_subplot(gs[1, :])

    # Rotate the spectra to match the metabolite map
    specs = np.rot90(spectra, axes=(0, 1))

    # Set default starting slice and coordinates if not provided
    if slice_idx is None:
        slice_idx = metab_map.shape[2] // 2
    if start_x is None:
        start_x = metab_map.shape[0] // 2
    if start_y is None:
        start_y = metab_map.shape[1] // 2

    # Plot the metabolite map
    im = ax_map.imshow(np.rot90(metab_map[:, :, slice_idx]), cmap='inferno', vmin=0)
    ax_map.set_title(f'{metab_name} Metabolite Map', color='white')
    ax_map.axis('off')

    # Add a colorbar to the left of the metabolite map
    cbar = fig.colorbar(im, ax=ax_map, orientation='vertical', pad=0.02, fraction=0.05)
    cbar.set_label('[mM / IU]', fontsize=14, weight='bold', color='white')
    cbar.ax.tick_params(axis='y', colors='white', labelsize=12)

    # Plot the initial MRSI spectra at the given starting coordinates
    spectrum_line, = ax_spec.plot(ppm_axis, specs[start_y, start_x, slice_idx].real, color='yellow')
    ax_spec.set_title(f'Spectrum Coordinates: ({start_x}, {start_y}, {slice_idx})', color='white')
    ax_spec.set_xlabel('ppm', color='white')
    ax_spec.set_ylabel('A.U.', color='white')
    ax_spec.set_xlim(5.5, -0.2)
    ax_spec.grid(True)

    # Add a green cross to the metabolite map to indicate the location of the MR spectrum
    cross, = ax_map.plot(start_x, start_y, 'o', markersize=16, markeredgewidth=4, markeredgecolor='green', markerfacecolor='none')

    # Set background color to black and label color to white
    fig.set_facecolor('black')
    ax_map.set_facecolor('black')
    ax_spec.set_facecolor('black')

    for ax in [ax_map, ax_spec]:
        ax.title.set_color('white')
        ax.xaxis.label.set_color('white')
        ax.yaxis.label.set_color('white')
        ax.tick_params(axis='x', colors='white')
        ax.tick_params(axis='y', colors='white')

    # Update slice, map, and spectrum when the slice index changes
    def update_slice(val):
        nonlocal slice_idx
        slice_idx = int(val)
        im.set_array(np.rot90(metab_map[:, :, slice_idx]))
        update_spectrum()

    # Update the spectrum at the current cross position
    def update_spectrum():
        x = int(cross.get_xdata()[0])
        y = int(cross.get_ydata()[0])
        new_spectrum = specs[y, x, slice_idx].real
        spectrum_line.set_ydata(new_spectrum)
        ax_spec.set_title(f'Spectrum at ({x}, {y}, {slice_idx})', color='white')

        # Adjust the y-axis limits dynamically
        min_y, max_y = new_spectrum.min(), new_spectrum.max()
        if min_y == max_y:
            range_y = 0.1 * abs(min_y) if min_y != 0 else 1
            ax_spec.set_ylim(min_y - range_y, max_y + range_y)
        else:
            ax_spec.set_ylim(min_y, max_y)

        fig.canvas.draw_idle()

    # Handle click events on the metabolite map
    def onclick(event):
        if event.inaxes == ax_map:
            x = int(event.xdata)
            y = int(event.ydata)
            cross.set_data([x], [y])
            update_spectrum()

    # Handle keypress events for slice navigation
    def onkeypress(event):
        nonlocal slice_idx
        if event.key == 'up' and slice_idx < metab_map.shape[2] - 1:
            slice_idx += 1
        elif event.key == 'down' and slice_idx > 0:
            slice_idx -= 1
        else:
            return
        slider.set_val(slice_idx)  # Sync slider with arrow keys
        update_slice(slice_idx)

    # Create the slider for slice navigation
    slider = Slider(ax_slider, 'Slice', 0, metab_map.shape[2] - 1, valinit=slice_idx, valstep=1, color='yellow')
    slider.on_changed(update_slice)

    # Connect events
    fig.canvas.mpl_connect('button_press_event', onclick)
    fig.canvas.mpl_connect('key_press_event', onkeypress)

    # Display the plot
    plt.show()

    return fig, (ax_map, ax_spec, ax_slider)







