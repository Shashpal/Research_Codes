#!/usr/bin/env python
# coding: utf-8

# In[1]:


from astropy.io import fits
import numpy as np

def subtract_fits(file1, file2, output_file):
    """
    Subtracts one FITS file from another and saves the result as a new FITS file.

    Parameters:
        file1 (str): Path to the first FITS file.
        file2 (str): Path to the second FITS file.
        output_file (str): Path to save the resulting FITS file.
    """
    # Open the first FITS file
    with fits.open(file1) as hdul1:
        data1 = hdul1[0].data  # Assuming image data is in the primary HDU
        header1 = hdul1[0].header  # Keep the header of the first file

    # Open the second FITS file
    with fits.open(file2) as hdul2:
        data2 = hdul2[0].data  # Assuming image data is in the primary HDU

    # Ensure both FITS files have the same dimensions
    if data1.shape != data2.shape:
        raise ValueError("The two FITS files have different dimensions and cannot be subtracted.")

    # Perform the subtraction
    result_data = data1 - data2

    # Create a new FITS HDU with the result
    hdu = fits.PrimaryHDU(data=result_data, header=header1)

    # Write the result to a new FITS file
    hdu.writeto(output_file, overwrite=True)
    print(f"Subtraction complete. Result saved to {output_file}")


# Example usage
file1 = "HE_0435/HE_0435_cropped.fits"  # Replace with the path to your first FITS file
file2 = "HE_0435/HE_0435_cropped_mask_residual_ring_for_SR.fits"  # Replace with the path to your second FITS file
output_file = "HE_0435/HE_0435_cropped_ring_SR.fits"  # Replace with the desired output file path

subtract_fits(file1, file2, output_file)

