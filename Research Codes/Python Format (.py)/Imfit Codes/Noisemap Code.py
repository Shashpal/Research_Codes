#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from astropy.io import fits
import os

def calculate_noisemap(fits_file, exposure_time, random_noise, output_file):
    """
    Calculate the photon noise and final noise map from a given FITS file.
    
    Parameters:
        fits_file (str): Path to the original FITS file.
        exposure_time (float): Exposure time in seconds.
        random_noise (float): Measured random noise of the image.
        output_file (str): Path to save the resulting noise map FITS file.
    """
    # Step 1: Load the original FITS file
    if not os.path.exists(fits_file):
        print(f"Error: File '{fits_file}' does not exist.")
        return
    
    with fits.open(fits_file) as hdul:
        original_data = hdul[1].data  # Assuming the FITS file has a single data extension
        header = hdul[1].header       # Keep the header for saving later
    
    # Step 2: Convert the FITS data to units of electrons by multiplying with exposure time
    data_in_electrons = original_data * exposure_time

    # Step 3: Select pixels with values larger than 3 * sigma of the random noise
    threshold = 3 * random_noise
    selected_pixels = np.where(data_in_electrons > threshold)

    # Step 4: Compute the photon noise for the selected pixels (sqrt of pixel value)
    photon_noise = np.zeros_like(data_in_electrons)  # Initialize an array with the same shape
    photon_noise[selected_pixels] = np.sqrt(data_in_electrons[selected_pixels])

    # Step 5: Convert back to units of electrons per second
    photon_noise_in_electrons_per_second = photon_noise / exposure_time

    # Step 6: Calculate the final noise map
    final_noise_map = np.sqrt(random_noise**2 + photon_noise_in_electrons_per_second**2)

    # Save the final noise map to a new FITS file
    hdu = fits.PrimaryHDU(data=final_noise_map, header=header)
    hdu.writeto(output_file, overwrite=True)
    print(f"Final noise map saved to '{output_file}'.")

# Example usage:
# Provide the path to your input FITS file, exposure time, random noise, and output file name
fits_file = "HE_0435/HE_0435_OG_Image.fits"  # Replace with your FITS file
exposure_time = 9940  # Replace with your exposure time (in seconds)
random_noise = 0.00738281  # Replace with your measured random noise (e.g., standard deviation of the background)
output_file = "HE_0435/noise_map_fixed.fits"  # Replace with your desired output file name

calculate_noisemap(fits_file, exposure_time, random_noise, output_file)

