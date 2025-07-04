import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider
import matplotlib.patches as patches
import numpy as np
import shutil

BASE_SIZE = 10

# Check if LaTeX is installed
latex_installed = shutil.which('latex') is not None
# print(f"Latex installed: {latex_installed}")

plt.rcParams.update({
    'text.usetex': latex_installed,
    'font.family': 'serif' if latex_installed else 'DejaVu Sans',
    'font.size': BASE_SIZE,
    'axes.titlesize': BASE_SIZE * 1.2,     
    'axes.labelsize': BASE_SIZE * 1.1,     
    'xtick.labelsize': BASE_SIZE * 0.9,    
    'ytick.labelsize': BASE_SIZE * 0.9,    
    'legend.fontsize': BASE_SIZE * 0.95,   
    'figure.titlesize': BASE_SIZE * 1.3,   
})


def get_slice(volume, slice_index=0, slice_axis=0):
    """
    Get a slice from a 3D volume along the specified axis.
    
    Args:
        volume (np.ndarray): 3D volume.
        slice_index (int): Index of the slice to extract.
        slice_axis (int): Axis along which to extract the slice (0, 1, or 2).
    
    Returns:
        np.ndarray: 2D slice of the volume.
    """
    if slice_axis == 0:
        return volume[slice_index, :, :]
    elif slice_axis == 1:
        return volume[:, slice_index, :]
    elif slice_axis == 2:
        return volume[:, :, slice_index]
    else:
        raise ValueError("slice_axis must be 0, 1, or 2.")
    

# Plotting
def plot_shim_map_process(label_map, base, boundary, final_map, slice_axis=2, slice_index=None, **kwargs):
    """
    Plot the shim map process.
    
    Args:
        label_map (np.ndarray): Label map.
        base (np.ndarray): Base shim map.
        boundary (np.ndarray): Boundary shim map.
        final_map (np.ndarray): Final shim map.
        slice_axis (int): Axis along which to extract the slice (0, 1, or 2).
        slice_index (int): Index of the slice to extract.
    """

    # Change orientation of label_map, base, boundary, and final_map to match orientation of the anatomical image
    label_map = label_map[::-1, ::-1, ::-1].copy()
    base = base[::-1, ::-1, ::-1].copy()
    boundary = boundary[::-1, ::-1, ::-1].copy()
    final_map = final_map[::-1, ::-1, ::-1].copy()

    if slice_index is None:
        slice_index = label_map.shape[slice_axis] // 2

    # Compute symmetric color range
    abs_max = max(
        abs(get_slice(base, slice_index, slice_axis)).max(),
        abs(get_slice(final_map, slice_index, slice_axis)).max()
    )
    vmin, vmax = -abs_max, abs_max

    # Create figure with GridSpec (3 rows, 3 columns)
    fig = plt.figure(figsize=(8, 8), dpi=300, constrained_layout=True)
    gs = gridspec.GridSpec(3, 3, height_ratios=[1, 1, 0.05], figure=fig)

    # Top row: 2 plots
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])

    # Bottom row: 2 plots
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    # Shared colorbar: occupies the same width as the 2x2 plot grid
    cax = fig.add_subplot(gs[2, 0:2])  # only the first two columns

    # Plotting
    ax1.imshow(get_slice(label_map, slice_index, slice_axis).T, cmap='grey')
    ax1.set_title('Label Map')

    ax2.imshow(get_slice(boundary, slice_index, slice_axis).T, cmap='grey')
    ax2.set_title('Boundary Map')

    im3 = ax3.imshow(get_slice(base, slice_index, slice_axis).T, cmap='bwr', vmin=vmin, vmax=vmax)
    ax3.set_title('Base Shim Map')

    im4 = ax4.imshow(get_slice(final_map, slice_index, slice_axis).T, cmap='bwr', vmin=vmin, vmax=vmax)
    ax4.set_title('Final Shim Map')

    # Clean look
    for ax in [ax1, ax2, ax3, ax4]:
        ax.axis('off')
        ax.set_aspect('equal')

    # Shared colorbar
    cbar = fig.colorbar(im4, cax=cax, orientation='horizontal')
    cbar.set_label('Shim Value (Hz)')

    if kwargs.get('show', False):
        plt.show()

    if 'save_path' in kwargs:
        plt.savefig(kwargs['save_path'], dpi=300, bbox_inches='tight')

def plot_spectra(ppm, specs, ppm_range=(0.5, 4.5), title=None, names=None, **kwargs):
    """
    Plot the spectrum.
    
    Args:
        spec (np.ndarray): Spectrum data.
        ppm (np.ndarray): PPM scale.
        ppm_range (tuple): Range of PPM values to display.
        title (str): Title of the plot.
    """
    plt.figure(figsize=(5, 5), dpi=300)
    for i, spec in enumerate(specs):
        if names is not None:
            plt.plot(ppm, spec, label=names[i])
        else:
            plt.plot(ppm, spec)
    plt.xlim(ppm_range[1], ppm_range[0])  # Inverted x-axis
    plt.grid()
    plt.legend()
    plt.xlabel('Chemical Shift (ppm)')
    plt.ylabel('Signal Intensity (a.u.)')
    plt.tight_layout()
    if title:
        plt.title(title)
    
    if kwargs.get('show', False):
        plt.show()

    if 'save_path' in kwargs:
        plt.savefig(kwargs['save_path'], dpi=300, bbox_inches='tight')

# Example inputs (replace these with your real data)
# img = np.random.rand(128, 128, 64)  # Replace with phantom.get_image_data()
# voi_coords = [30, 60, 40, 70, 10, 40]  # Replace with your VOI selection

def plot_voi_scroll(img, voi_coords):
    x_min, x_max, y_min, y_max, z_min, z_max = voi_coords

    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    plt.subplots_adjust(bottom=0.25)

    # Initial slice indices
    slice_z = (z_min + z_max) // 2
    slice_y = (y_min + y_max) // 2
    slice_x = (x_min + x_max) // 2

    # Display axial (xy at fixed z)
    im0 = axs[0].imshow(img[:, :, slice_z].T, origin='lower', cmap='gray')
    rect0 = patches.Rectangle((x_min, y_min), x_max - x_min, y_max - y_min,
                              linewidth=1.5, edgecolor='r', facecolor='none')
    axs[0].add_patch(rect0)
    axs[0].set_title(f'Axial (z={slice_z})')

    # Display coronal (xz at fixed y)
    im1 = axs[1].imshow(img[:, slice_y, :].T, origin='lower', cmap='gray')
    rect1 = patches.Rectangle((x_min, z_min), x_max - x_min, z_max - z_min,
                              linewidth=1.5, edgecolor='r', facecolor='none')
    axs[1].add_patch(rect1)
    axs[1].set_title(f'Coronal (y={slice_y})')

    # Display sagittal (yz at fixed x)
    im2 = axs[2].imshow(img[slice_x, :, :].T, origin='lower', cmap='gray')
    rect2 = patches.Rectangle((y_min, z_min), y_max - y_min, z_max - z_min,
                              linewidth=1.5, edgecolor='r', facecolor='none')
    axs[2].add_patch(rect2)
    axs[2].set_title(f'Sagittal (x={slice_x})')

    for ax in axs:
        ax.set_xticks([])
        ax.set_yticks([])

    # Add sliders for z, y, x
    axcolor = 'lightgoldenrodyellow'
    ax_z = plt.axes([0.15, 0.15, 0.65, 0.03], facecolor=axcolor)
    ax_y = plt.axes([0.15, 0.10, 0.65, 0.03], facecolor=axcolor)
    ax_x = plt.axes([0.15, 0.05, 0.65, 0.03], facecolor=axcolor)

    slider_z = Slider(ax_z, 'Z', 0, img.shape[2] - 1, valinit=slice_z, valfmt='%d')
    slider_y = Slider(ax_y, 'Y', 0, img.shape[1] - 1, valinit=slice_y, valfmt='%d')
    slider_x = Slider(ax_x, 'X', 0, img.shape[0] - 1, valinit=slice_x, valfmt='%d')

    def update(val):
        z = int(slider_z.val)
        y = int(slider_y.val)
        x = int(slider_x.val)

        im0.set_data(img[:, :, z].T)
        axs[0].set_title(f'Axial (z={z})')

        im1.set_data(img[:, y, :].T)
        axs[1].set_title(f'Coronal (y={y})')

        im2.set_data(img[x, :, :].T)
        axs[2].set_title(f'Sagittal (x={x})')

        fig.canvas.draw_idle()

    slider_z.on_changed(update)
    slider_y.on_changed(update)
    slider_x.on_changed(update)

    plt.show()

    import numpy as np
import matplotlib.pyplot as plt

def plot_avg_spectra_with_std_and_diff(specs_real_norm, specs_sim_norm, ppm_axis, ppm_range, label_real="In-vivo", label_sim="Simulated"):
    """
    Plots separated average spectra with shaded std (real and simulated), and their difference.

    Parameters:
    - specs_real_norm: np.ndarray, shape (n_real, n_points), complex or real
    - specs_sim_norm: np.ndarray, shape (n_sim, n_points), complex or real
    - ppm_axis: np.ndarray, shape (n_points,)
    - ppm_range: tuple (min_ppm, max_ppm)
    - label_real: str, label for real spectra
    - label_sim: str, label for simulated spectra
    """

    def find_ppm_indices(ppm, start, end):
        return np.where((ppm >= start) & (ppm <= end))[0]

    def extract_real_range(spectra, ppm_inds):
        return np.real(spectra[:, ppm_inds])

    # Find indices and select range
    pca_inds = find_ppm_indices(ppm_axis, *ppm_range)
    ppm_selected = ppm_axis[pca_inds]

    # Extract real parts in range
    real_spectra = extract_real_range(specs_real_norm, pca_inds)
    sim_spectra = extract_real_range(specs_sim_norm, pca_inds)

    # Compute averages and stds
    real_mean, real_std = np.mean(real_spectra, axis=0), np.std(real_spectra, axis=0)
    sim_mean, sim_std = np.mean(sim_spectra, axis=0), np.std(sim_spectra, axis=0)
    diff = real_mean - sim_mean
    diff_std = np.sqrt(real_std**2 + sim_std**2)  # Combined std for difference

    # Plot layout: 3 vertically stacked subplots
    fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True, gridspec_kw={'height_ratios': [1, 1, 0.6]})

    # Plot In-vivo
    axs[0].plot(ppm_selected, real_mean, color='tab:blue', label=f'{label_real} average')
    axs[0].fill_between(ppm_selected, real_mean - real_std, real_mean + real_std, color='tab:blue', alpha=0.3)
    axs[0].set_ylabel('Signal (a.u.)')
    axs[0].legend()
    axs[0].set_title('Average MRS Spectra')

    # Plot Simulated
    axs[1].plot(ppm_selected, sim_mean, color='tab:red', label=f'{label_sim} average')
    axs[1].fill_between(ppm_selected, sim_mean - sim_std, sim_mean + sim_std, color='tab:red', alpha=0.3)
    axs[1].set_ylabel('Signal (a.u.)')
    axs[1].legend()

    # Plot Difference
    axs[2].plot(ppm_selected, diff, color='tab:green', label='Difference')
    axs[2].fill_between(ppm_selected, diff - diff_std, diff + diff_std, color='tab:green', alpha=0.3)
    axs[2].axhline(0, color='black', linestyle='--', linewidth=1)
    axs[2].set_ylabel('Difference')
    axs[2].set_xlabel('ppm')
    axs[2].set_title(f'Difference: {label_real} - {label_sim}')

    # Aesthetics
    for ax in axs:
        ax.grid(True)
    axs[-1].invert_xaxis()  # ppm axis decreases left to right
    plt.tight_layout()
    plt.show()
