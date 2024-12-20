#!/bin/bash

# Initialize an array to track failed packages
failed_packages=()

# Function to install a package and track failures
install_with_conda() {
    package=$1
    echo "Installing: $package"
    if ! conda install -y $package; then
        echo "Failed to install: $package"
        failed_packages+=("$package")
    fi
}

# Function to install a package from specific channels
install_with_conda_channels() {
    package=$1
    channels=$2
    echo "Installing: $package from $channels"
    if ! conda install -y $channels $package; then
        echo "Failed to install: $package"
        failed_packages+=("$package")
    fi
}

# Update conda
echo "Updating conda..."
if ! conda update -n base -c defaults conda -y; then
    echo "Warning: Conda update failed, continuing with installation..."
fi

# Install packages
install_with_conda "ipykernel"
install_with_conda_channels "fsl_mrs" "-c conda-forge -c https://fsl.fmrib.ox.ac.uk/fsldownloads/fslconda/public/"
install_with_conda_channels "torchio" "-c conda-forge"
install_with_conda_channels "ipywidgets" "-c conda-forge"
install_with_conda "openpyxl"
install_with_conda_channels "ipympl" "-c conda-forge"
install_with_conda_channels "jupyterlab" "-c conda-forge"

# Print summary
echo "Installation Summary:"
if [ ${#failed_packages[@]} -eq 0 ]; then
    echo "All packages installed successfully!"
else
    echo "The following packages failed to install:"
    for pkg in "${failed_packages[@]}"; do
        echo "- $pkg"
    done
fi
