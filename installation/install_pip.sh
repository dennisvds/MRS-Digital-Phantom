#!/bin/bash

# Initialize an array to track failed packages
failed_packages=()

# Function to install a package and track failures
install_with_pip() {
    package=$1
    echo "Installing: $package"
    if ! pip install $package; then
        echo "Failed to install: $package"
        failed_packages+=("$package")
    fi
}

# Function to install FSL-specific package (for now, we'll keep it out of pip due to channel-specific requirements)
install_fsl_mrs() {
    package=$1
    echo "Installing: $package from specific FSL channel"
    # Assuming FSL is available via pip, but if not, we'll need a different solution.
    if ! pip install $package; then
        echo "Failed to install: $package"
        failed_packages+=("$package")
    fi
}

# Update pip
echo "Updating pip..."
if ! pip install --upgrade pip; then
    echo "Warning: Pip update failed, continuing with installation..."
fi

# Install packages
install_with_pip "ipykernel"
install_fsl_mrs "fsl-mrs"  # If FSL-MRS is available on PyPI; otherwise, this can be a custom installation process
install_with_pip "torchio"
install_with_pip "ipywidgets"
install_with_pip "openpyxl"
install_with_pip "ipympl"
install_with_pip "jupyterlab"

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
