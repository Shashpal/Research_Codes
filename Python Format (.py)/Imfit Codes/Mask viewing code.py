#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
from astropy.io import fits

def apply_mask(original_fits, mask_fits, output_fits):
    """
    Apply a mask to an original FITS file and create a new FITS file
    showing the leftover light after masking.

    Parameters:
    - original_fits: str, path to the original FITS file containing all the light.
    - mask_fits: str, path to the mask FITS file (1 for masked regions, 0 for unmasked).
    - output_fits: str, path to save the output FITS file.

    Returns:
    - None. The result is saved as a new FITS file.
    """
    # Load the original FITS file
    with fits.open(original_fits) as hdul_original:
        original_data = hdul_original[0].data
        original_header = hdul_original[0].header

    # Load the mask FITS file
    with fits.open(mask_fits) as hdul_mask:
        mask_data = hdul_mask[0].data

    # Sanity check: Ensure the dimensions of the two files match
    if original_data.shape != mask_data.shape:
        raise ValueError("The dimensions of the original FITS file and the mask FITS file do not match.")

    # Apply the mask: Set regions where mask_data == 1 to 0 in the original data
    masked_data = np.where(mask_data == 1, 0, original_data)

    # Save the masked data to a new FITS file
    hdu = fits.PrimaryHDU(masked_data, header=original_header)
    hdu.writeto(output_fits, overwrite=True)

    print(f"Masked FITS file created: {output_fits}")


# Example usage
if __name__ == "__main__":
    # Path to the original FITS file
    original_fits = "HE_0435/HE_0435_cropped_residual_ring.fits"  # Replace with your original FITS file path

    # Path to the mask FITS file
    mask_fits = "HE_0435/HE_0435_cropped_mask_residual_ring.fits"  # Replace with your mask FITS file path

    # Path to save the output FITS file
    output_fits = "HE_0435/HE_0435_cropped_mask_residual_ring_for_SR.fits"  # Replace with desired output path

    # Run the function
    apply_mask(original_fits, mask_fits, output_fits)

