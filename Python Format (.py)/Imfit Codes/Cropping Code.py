#!/usr/bin/env python
# coding: utf-8

# In[3]:


# Code to crop the fits file to a smaller fits file 
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from regions import Regions

def crop_fits(fits_file, region_file, output_fits_file, hdu_index=0):
    """
    Crop a FITS file based on a region file. The output FITS file will only contain
    data within the specified region, with all pixels outside the region set to NaN.

    Parameters:
    - fits_file: str, path to the input FITS file.
    - region_file: str, path to the DS9 region file.
    - output_fits_file: str, path to save the cropped FITS file.
    - hdu_index: int, index of the HDU to use from the FITS file (default is 0).

    Returns:
    - None. The cropped FITS file is saved to output_fits_file.
    """

    # Load the FITS file
    with fits.open(fits_file) as hdul:
        # Check if the provided HDU index is valid
        if hdu_index >= len(hdul):
            raise ValueError(f"HDU index {hdu_index} is out of range. FITS file has {len(hdul)} HDUs.")

        # Select the specified HDU
        hdu = hdul[hdu_index]
        if hdu.data is None:
            raise ValueError(f"The selected HDU[{hdu_index}] does not contain image data.")

        # Extract data and header from the selected HDU
        data = hdu.data
        header = hdu.header

        # Load the WCS from the header
        wcs = WCS(header)

    # Load the region file
    regions = Regions.read(region_file, format="ds9")

    # Check if the region file contains at least one region
    if len(regions) == 0:
        raise ValueError("No regions found in the provided region file.")

    # Use the first region in the file (assuming only one region for cropping)
    region = regions[0]

    # Convert the region to pixel coordinates (if needed)
    if hasattr(region, "to_pixel"):
        pixel_region = region.to_pixel(wcs)
    else:
        pixel_region = region

    # Create a mask for the region
    region_mask = pixel_region.to_mask(mode="center")

    # Crop the data using the region mask
    cropped_data = region_mask.cutout(data)

    # Check if the cutout is valid
    if cropped_data is None:
        raise ValueError("The region is completely outside the bounds of the FITS file.")

    # Apply the mask to the cropped data
    mask = region_mask.to_image(region_mask.shape)
    if mask is not None:  # Ensure the mask is valid
        cropped_data[mask == 0] = np.nan  # Set pixels outside the region to NaN

    # Update the header to reflect the cropped region
    bb = pixel_region.bounding_box
    new_header = header.copy()
    new_header['CRPIX1'] -= bb.ixmin
    new_header['CRPIX2'] -= bb.iymin

    # Save the cropped data as a new FITS file
    hdu = fits.PrimaryHDU(cropped_data, header=new_header)
    hdu.writeto(output_fits_file, overwrite=True)

    print(f"Cropped FITS file created: {output_fits_file}")


# Example usage
if __name__ == "__main__":
    # Input FITS file
    fits_file = "HE_0435/noise_map.fits"  # Input original FITS file

    # Input region file
    region_file = "HE_0435/HE_0435 Crop.reg"  # Replace with your region file path

    # Output cropped FITS file
    output_fits_file = "HE_0435/noise_map_cropped.fits"  # Replace with desired output path

    # Specify the HDU index you want to use
    hdu_index = 0  # Replace with the desired HDU index (e.g., 0, 1, 2, etc.)

    # Run the function
    crop_fits(fits_file, region_file, output_fits_file, hdu_index=hdu_index)

