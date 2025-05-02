#!/usr/bin/env python
# coding: utf-8

# In[1]:


from astropy.io import fits
fits.info("HE_0435/HE_0435_OG_Image.fits")

# i.e Dimension of (1391, 1268) refers to 1391 pixels wide (columns) and 1268 pixels tall (rows)


# In[5]:


# Code to mask pixels in the region to have a pixel value of 1 
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from regions import Regions, CircleAnnulusPixelRegion, PixCoord

def create_mask(fits_file, region_file, output_mask_file, hdu_index=0):
    """
    Create a masking FITS file based on a region file. The mask will have
    pixel values of 1 inside the regions and 0 outside.

    Parameters:
    - fits_file: str, path to the input FITS file.
    - region_file: str, path to the DS9 region file.
    - output_mask_file: str, path to save the output masking FITS file.
    - hdu_index: int, index of the HDU to use from the FITS file (default is 0).

    Returns:
    - None. The masking FITS file is saved to output_mask_file.
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

    # Create an empty mask with the same shape as the FITS file
    mask = np.zeros(data.shape, dtype=np.int16)

    # Load the region file
    regions = Regions.read(region_file, format="ds9")

    # Process each region
    for region in regions:
        # Check if the region is in pixel or sky coordinates
        if hasattr(region, "to_pixel"):
            # Convert sky regions to pixel regions
            pixel_region = region.to_pixel(wcs)
        else:
            # Region is already in pixel coordinates
            pixel_region = region

        # Create a mask for the current region
        region_mask = pixel_region.to_mask(mode="center")
        mask_data = region_mask.to_image(data.shape)

        # Add the region mask to the overall mask
        if mask_data is not None:  # Ensure the region is within the image bounds
            mask += (mask_data > 0).astype(np.int16)

    # Ensure the mask only contains 0 and 1 values
    mask = np.clip(mask, 0, 1)

    # Save the mask as a new FITS file
    hdu = fits.PrimaryHDU(mask, header=header)
    hdu.writeto(output_mask_file, overwrite=True)

    print(f"Masking FITS file created: {output_mask_file}")


# Example usage
if __name__ == "__main__":
    # Input FITS file
    fits_file = "HE_0435/HE_0435_cropped_ring_SR.fits"  # Input original fits file 

    # Input region file
    region_file = "HE_0435/HE_0435 Mask_Halo_SR.reg"  # Replace with your region file path

    # Output mask FITS file
    output_mask_file = "HE_0435/HE_0435_cropped_mask_ring_Halo_SR.fits"  # Replace with desired output path

    # Specify the HDU index you want to use
    hdu_index = 0  # Replace with the desired HDU index (e.g., 0, 1, 2, etc.)

    # Run the function
    create_mask(fits_file, region_file, output_mask_file, hdu_index=hdu_index)


# In[ ]:




