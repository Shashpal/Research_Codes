#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Code that completely removes pairs based on the observed flux ratio. Used for when the parity 
# normalisation isn't 1.0 flux ratio
import numpy as np
from astropy.io import fits

# Filenames for 1 input constraint file and multiple fits and output file pairs 
Input_file = "RXJ0911/SPC/obs_point.dat"
Output_files = [
    "RXJ0911/SPC/outSP_point.dat",
    "RXJ0911/SPC/outSPG_point.dat",
    "RXJ0911/SPC/outSPR_point.dat",
    "RXJ0911/SPC/outSPGR_point.dat",
    "RXJ0911/SPFC/outSPF_point.dat",
    "RXJ0911/SPFC/outSPFG_point.dat",
    "RXJ0911/SPFC/outSPFR_point.dat",
    "RXJ0911/SPFC/outSPFGR_point.dat",
    "RXJ0911/PPC/outPP_point.dat",
    "RXJ0911/PPC/outPPR_point.dat",
    "RXJ0911/PPC/outPPG_point.dat",
    "RXJ0911/PPC/outPPGR_point.dat",
    "RXJ0911/PPFC/outPPF_point.dat",
    "RXJ0911/PPFC/outPPFR_point.dat",
    "RXJ0911/PPFC/outPPFG_point.dat",
    "RXJ0911/PPFC/outPPFGR_point.dat",
    "RXJ0911/NPC/outNP_point.dat",
    "RXJ0911/NPC/outNPR_point.dat",
    "RXJ0911/NPC/outNPG_point.dat",
    "RXJ0911/NPC/outNPGR_point.dat",
    "RXJ0911/NPFC/outNPF_point.dat",
    "RXJ0911/NPFC/outNPFR_point.dat",
    "RXJ0911/NPFC/outNPFG_point.dat",
    "RXJ0911/NPFC/outNPFGR_point.dat"
]
Fits_files = [
    "RXJ0911/SPC/outSP_lens.fits",
    "RXJ0911/SPC/outSPG_lens.fits",
    "RXJ0911/SPC/outSPR_lens.fits",
    "RXJ0911/SPC/outSPGR_lens.fits",
    "RXJ0911/SPFC/outSPF_lens.fits",
    "RXJ0911/SPFC/outSPFG_lens.fits",
    "RXJ0911/SPFC/outSPFR_lens.fits",
    "RXJ0911/SPFC/outSPFGR_lens.fits",
    "RXJ0911/PPC/outPP_lens.fits",
    "RXJ0911/PPC/outPPR_lens.fits",
    "RXJ0911/PPC/outPPG_lens.fits",
    "RXJ0911/PPC/outPPGR_lens.fits",
    "RXJ0911/PPFC/outPPF_lens.fits",
    "RXJ0911/PPFC/outPPFR_lens.fits",
    "RXJ0911/PPFC/outPPFG_lens.fits",
    "RXJ0911/PPFC/outPPFGR_lens.fits",
    "RXJ0911/NPC/outNP_lens.fits",
    "RXJ0911/NPC/outNPR_lens.fits",
    "RXJ0911/NPC/outNPG_lens.fits",
    "RXJ0911/NPC/outNPGR_lens.fits",
    "RXJ0911/NPFC/outNPF_lens.fits",
    "RXJ0911/NPFC/outNPFR_lens.fits",
    "RXJ0911/NPFC/outNPFG_lens.fits",
    "RXJ0911/NPFC/outNPFGR_lens.fits"
]
Output_text_file = "RXJ0911_Parity_4thB_test.txt"

n = 0.01  # Scaling factor for the coordinates-to-pixel conversion

# Function to calculate distance between two points
def calculate_distance(x1, y1, x2, y2):
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

# Read observed image positions (x1, y1) and flux ratios from obs_point.dat
observed_positions = []
observed_flux_ratios = [0.56, 1, 0.24, 0.53]  # Flux ratios at observed positions
with open(Input_file, 'r') as obs_file:
    lines = obs_file.readlines()[1:]  # Skip the first line
    for line in lines:
        columns = line.split()
        if len(columns) >= 2:
            x1, y1 = float(columns[0]), float(columns[1])
            observed_positions.append((x1, y1))

# Initialize statistical counters
negative_parity_pred_yes = 0
negative_parity_pred_total = 0
positive_parity_pred_yes = 0
positive_parity_pred_total = 0

negative_parity_obs_yes = 0
negative_parity_obs_total = 0
positive_parity_obs_yes = 0
positive_parity_obs_total = 0

yes_pred_count = 0
no_pred_count = 0
yes_obs_count = 0
no_obs_count = 0

# Open the output file in write mode
with open(Output_text_file, 'w') as output_file:

    # Loop over each model (Output_file and Fits_file pairs)
    for output_filename, fits_filename in zip(Output_files, Fits_files):
        # Read predicted image positions (x2, y2) and magnifications (M2) from output file
        predicted_positions = []
        predicted_magnifications = []
        original_magnifications = []  # Store the non-absolute magnifications
        with open(output_filename, 'r') as out_file:
            lines = out_file.readlines()[1:]  # Skip the first line
            for line in lines:
                columns = line.split()
                if len(columns) >= 4 and not columns[0].startswith('#'):
                    try:
                        x2, y2 = float(columns[0]), float(columns[1])
                        magnification = float(columns[2])  # Store original magnification (signed)
                        predicted_positions.append((x2, y2))
                        predicted_magnifications.append(abs(magnification))  # Store absolute magnification
                        original_magnifications.append(magnification)
                    except ValueError:
                        pass

        # Pair observed and predicted positions based on minimum distance
        paired_positions = []
        paired_magnifications = []
        paired_original_magnifications = []
        for obs_pos in observed_positions:
            min_distance = float('inf')
            closest_pred_pos = None
            closest_magnification = None
            closest_original_magnification = None
            for pred_pos, magnification, orig_magnification in zip(predicted_positions, predicted_magnifications, original_magnifications):
                distance = calculate_distance(obs_pos[0], obs_pos[1], pred_pos[0], pred_pos[1])
                if distance < min_distance:
                    min_distance = distance
                    closest_pred_pos = pred_pos
                    closest_magnification = magnification
                    closest_original_magnification = orig_magnification
            paired_positions.append((obs_pos, closest_pred_pos))
            paired_magnifications.append(closest_magnification)
            paired_original_magnifications.append(closest_original_magnification)

        # Open FITS file to extract image magnifications
        with fits.open(fits_filename) as hdul:
            Mag = hdul[0].data[6]  # Assuming the magnification data is in the 6th HDU

        # Convert observed positions to pixel coordinates and get magnifications from the FITS file
        fits_magnifications = []
        image_size = Mag.shape
        center_x = image_size[1] // 2
        center_y = image_size[0] // 2

        # Convert observed positions to pixel coordinates and retrieve magnifications
        pixel_coords = []
        for obs_pos in observed_positions:
            pixel_x = int((obs_pos[0] + center_x * n) / n)
            pixel_y = int((obs_pos[1] + center_y * n) / n)
            pixel_coords.append((pixel_x, pixel_y))

            try:
                fits_magnification = Mag[pixel_y, pixel_x]  # Get magnification from FITS file
            except IndexError:
                fits_magnification = np.nan  # Handle out-of-bounds access
            fits_magnifications.append(fits_magnification)

        # Calculate the predicted flux ratios at the predicted positions
        green_image_index = 2  # Image order used as a flux ratio denominator
        predicted_flux_ratios = [magnification / predicted_magnifications[green_image_index] for magnification in predicted_magnifications]

        # Calculate the predicted flux ratios at the observed positions using the FITS file
        abs_values = [1 / abs(value) if value != 0 else np.nan for value in fits_magnifications]
        predicted_flux_ratios_observed = [abs_values[i] / abs_values[green_image_index] for i in range(len(abs_values))]

        # Also divide the observed flux ratios
        observed_flux_ratios_adjusted = [observed_flux_ratios[i] / observed_flux_ratios[green_image_index] for i in range(len(observed_flux_ratios))]

        # Print results for the current model
        output_file.write(f"\nModel: {output_filename}\n")
        print(f"\nModel: {output_filename}\n")

        for i, (obs_pos, pred_pos) in enumerate(paired_positions): 
            # Exclude pairs with an observed flux ratio of 1.0
            if observed_flux_ratios[i] == 0.24: 
                continue

            pred_sign = "(pred_negative)" if paired_original_magnifications[i] < 0 else "(pred_positive)"
            obs_sign = "(obs_negative)" if fits_magnifications[i] < 0 else "(obs_positive)"

            # Check if the predicted flux ratio is higher or lower
            pred_phrase = "(over_pred)" if predicted_flux_ratios[i] > observed_flux_ratios_adjusted[i] else "(under_pred)"
            obs_phrase = "(over_obs)" if predicted_flux_ratios_observed[i] > observed_flux_ratios_adjusted[i] else "(under_obs)"

            # Add YES/NO tags
            pred_tag = "[YES_pred]" if (pred_sign == "(pred_negative)" and pred_phrase == "(over_pred)") or                                        (pred_sign == "(pred_positive)" and pred_phrase == "(under_pred)") else "[NO_pred]"
            obs_tag = "[YES_obs]" if (obs_sign == "(obs_negative)" and obs_phrase == "(over_obs)") or                                    (obs_sign == "(obs_positive)" and obs_phrase == "(under_obs)") else "[NO_obs]"

            result_str = (
                f"Pair {i+1} {pred_sign} {obs_sign} {pred_phrase} {obs_phrase} {pred_tag} {obs_tag}:\n"
                f"Observed Position: {obs_pos}\n"
                f"Closest Predicted Position: {pred_pos}\n"
                f"Observed Flux Ratio: {observed_flux_ratios_adjusted[i]:.2f}\n"
                f"Predicted Magnification: {paired_original_magnifications[i]:.2f}\n"
                f"Magnification from FITS file: {fits_magnifications[i]:.2f}\n"
                f"Predicted Flux Ratio at Predicted Position: {predicted_flux_ratios[i]:.2f}\n"
                f"Predicted Flux Ratio at Observed Position (FITS): {predicted_flux_ratios_observed[i]:.2f}\n\n"
            )

            output_file.write(result_str)
            print(result_str)

            # Update counters for statistical analysis
            if pred_sign == "(pred_negative)":
                negative_parity_pred_total += 1
                if pred_tag == "[YES_pred]":
                    negative_parity_pred_yes += 1
            elif pred_sign == "(pred_positive)":
                positive_parity_pred_total += 1
                if pred_tag == "[YES_pred]":
                    positive_parity_pred_yes += 1

            if obs_sign == "(obs_negative)":
                negative_parity_obs_total += 1
                if obs_tag == "[YES_obs]":
                    negative_parity_obs_yes += 1
            elif obs_sign == "(obs_positive)":
                positive_parity_obs_total += 1
                if obs_tag == "[YES_obs]":
                    positive_parity_obs_yes += 1

            if pred_tag == "[YES_pred]":
                yes_pred_count += 1
            elif pred_tag == "[NO_pred]":
                no_pred_count += 1

            if obs_tag == "[YES_obs]":
                yes_obs_count += 1
            elif obs_tag == "[NO_obs]":
                no_obs_count += 1

# Calculate statistical analysis results
def calculate_percentage(numerator, denominator):
    return (numerator / denominator * 100) if denominator > 0 else 0

# Predicted positions
negative_parity_pred_benefit = calculate_percentage(negative_parity_pred_yes, negative_parity_pred_total)
positive_parity_pred_benefit = calculate_percentage(positive_parity_pred_yes, positive_parity_pred_total)

# Observed positions
negative_parity_obs_benefit = calculate_percentage(negative_parity_obs_yes, negative_parity_obs_total)
positive_parity_obs_benefit = calculate_percentage(positive_parity_obs_yes, positive_parity_obs_total)

# Overall benefit
overall_pred_benefit = calculate_percentage(yes_pred_count, yes_pred_count + no_pred_count)
overall_obs_benefit = calculate_percentage(yes_obs_count, yes_obs_count + no_obs_count)

# Output the statistical results in the requested format
with open(Output_text_file, 'a') as output_file:
    results = (
        "\nStatistical Analysis:\n"
        f"% of negative parity that will benefit from adding subhaloes (predicted positions): "
        f"{negative_parity_pred_yes}/{negative_parity_pred_total} x 100% = {negative_parity_pred_benefit:.2f}%\n"
        f"% of positive parity that will benefit from adding subhaloes (predicted positions): "
        f"{positive_parity_pred_yes}/{positive_parity_pred_total} x 100% = {positive_parity_pred_benefit:.2f}%\n"
        f"% of negative parity that will benefit from adding subhaloes (observed positions): "
        f"{negative_parity_obs_yes}/{negative_parity_obs_total} x 100% = {negative_parity_obs_benefit:.2f}%\n"
        f"% of positive parity that will benefit from adding subhaloes (observed positions): "
        f"{positive_parity_obs_yes}/{positive_parity_obs_total} x 100% = {positive_parity_obs_benefit:.2f}%\n"
        f"% of images benefited from adding subhaloes (predicted positions): "
        f"{yes_pred_count}/{yes_pred_count + no_pred_count} x 100% = {overall_pred_benefit:.2f}%\n"
        f"% of images benefited from adding subhaloes (observed positions): "
        f"{yes_obs_count}/{yes_obs_count + no_obs_count} x 100% = {overall_obs_benefit:.2f}%\n"
    )

    output_file.write(results)
    print(results)


# In[ ]:


# Code that properly ignores pairs that don't have [YES_pred] and so on. It still displays all 4 pairs
# so it is suitable for PSJ1606 which has two 1.0 flux ratio models
import numpy as np
from astropy.io import fits

# Filenames for 1 input constraint file and multiple fits and output file pairs 
Input_file = "RXJ0911/SPC/obs_point.dat"
Output_files = [
    "RXJ0911/SPC/outSP_point.dat",
    "RXJ0911/SPC/outSPG_point.dat",
    "RXJ0911/SPC/outSPR_point.dat",
    "RXJ0911/SPC/outSPGR_point.dat",
    "RXJ0911/SPFC/outSPF_point.dat",
    "RXJ0911/SPFC/outSPFG_point.dat",
    "RXJ0911/SPFC/outSPFR_point.dat",
    "RXJ0911/SPFC/outSPFGR_point.dat",
    "RXJ0911/PPC/outPP_point.dat",
    "RXJ0911/PPC/outPPR_point.dat",
    "RXJ0911/PPC/outPPG_point.dat",
    "RXJ0911/PPC/outPPGR_point.dat",
    "RXJ0911/PPFC/outPPF_point.dat",
    "RXJ0911/PPFC/outPPFR_point.dat",
    "RXJ0911/PPFC/outPPFG_point.dat",
    "RXJ0911/PPFC/outPPFGR_point.dat",
    "RXJ0911/NPC/outNP_point.dat",
    "RXJ0911/NPC/outNPR_point.dat",
    "RXJ0911/NPC/outNPG_point.dat",
    "RXJ0911/NPC/outNPGR_point.dat",
    "RXJ0911/NPFC/outNPF_point.dat",
    "RXJ0911/NPFC/outNPFR_point.dat",
    "RXJ0911/NPFC/outNPFG_point.dat",
    "RXJ0911/NPFC/outNPFGR_point.dat"
]
Fits_files = [
    "RXJ0911/SPC/outSP_lens.fits",
    "RXJ0911/SPC/outSPG_lens.fits",
    "RXJ0911/SPC/outSPR_lens.fits",
    "RXJ0911/SPC/outSPGR_lens.fits",
    "RXJ0911/SPFC/outSPF_lens.fits",
    "RXJ0911/SPFC/outSPFG_lens.fits",
    "RXJ0911/SPFC/outSPFR_lens.fits",
    "RXJ0911/SPFC/outSPFGR_lens.fits",
    "RXJ0911/PPC/outPP_lens.fits",
    "RXJ0911/PPC/outPPR_lens.fits",
    "RXJ0911/PPC/outPPG_lens.fits",
    "RXJ0911/PPC/outPPGR_lens.fits",
    "RXJ0911/PPFC/outPPF_lens.fits",
    "RXJ0911/PPFC/outPPFR_lens.fits",
    "RXJ0911/PPFC/outPPFG_lens.fits",
    "RXJ0911/PPFC/outPPFGR_lens.fits",
    "RXJ0911/NPC/outNP_lens.fits",
    "RXJ0911/NPC/outNPR_lens.fits",
    "RXJ0911/NPC/outNPG_lens.fits",
    "RXJ0911/NPC/outNPGR_lens.fits",
    "RXJ0911/NPFC/outNPF_lens.fits",
    "RXJ0911/NPFC/outNPFR_lens.fits",
    "RXJ0911/NPFC/outNPFG_lens.fits",
    "RXJ0911/NPFC/outNPFGR_lens.fits"
]
Output_text_file = "RXJ0911_Parity_1stB_test.txt"

n = 0.01

# Function to calculate distance between two points
def calculate_distance(x1, y1, x2, y2):
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

# Read observed image positions (x1, y1) from obs_point.dat
observed_positions = []
observed_flux_ratios = [0.56, 1, 0.24, 0.53]
with open(Input_file, 'r') as obs_file:
    lines = obs_file.readlines()[1:]
    for line in lines:
        columns = line.split()
        if len(columns) >= 2:
            x1, y1 = float(columns[0]), float(columns[1])
            observed_positions.append((x1, y1))

# Initialize statistical counters
negative_parity_pred_yes = 0
negative_parity_pred_total = 0
positive_parity_pred_yes = 0
positive_parity_pred_total = 0

negative_parity_obs_yes = 0
negative_parity_obs_total = 0
positive_parity_obs_yes = 0
positive_parity_obs_total = 0

yes_pred_count = 0
no_pred_count = 0
yes_obs_count = 0
no_obs_count = 0

# Open the output file in write mode
with open(Output_text_file, 'w') as output_file:

    # Loop over each model (Output_file and Fits_file pairs)
    for output_filename, fits_filename in zip(Output_files, Fits_files):
        predicted_positions = []
        predicted_magnifications = []
        original_magnifications = []
        with open(output_filename, 'r') as out_file:
            lines = out_file.readlines()[1:]
            for line in lines:
                columns = line.split()
                if len(columns) >= 4 and not columns[0].startswith('#'):
                    try:
                        x2, y2 = float(columns[0]), float(columns[1])
                        magnification = float(columns[2])
                        predicted_positions.append((x2, y2))
                        predicted_magnifications.append(abs(magnification))
                        original_magnifications.append(magnification)
                    except ValueError:
                        pass

        # Pair observed and predicted positions based on minimum distance
        paired_positions = []
        paired_magnifications = []
        paired_original_magnifications = []
        for obs_pos in observed_positions:
            min_distance = float('inf')
            closest_pred_pos = None
            closest_magnification = None
            closest_original_magnification = None
            for pred_pos, magnification, orig_magnification in zip(predicted_positions, predicted_magnifications, original_magnifications):
                distance = calculate_distance(obs_pos[0], obs_pos[1], pred_pos[0], pred_pos[1])
                if distance < min_distance:
                    min_distance = distance
                    closest_pred_pos = pred_pos
                    closest_magnification = magnification
                    closest_original_magnification = orig_magnification
            paired_positions.append((obs_pos, closest_pred_pos))
            paired_magnifications.append(closest_magnification)
            paired_original_magnifications.append(closest_original_magnification)

        # Open FITS file to extract image magnifications
        with fits.open(fits_filename) as hdul:
            Mag = hdul[0].data[6]

        # Convert observed positions to pixel coordinates and get magnifications from the FITS file
        fits_magnifications = []
        image_size = Mag.shape
        center_x = image_size[1] // 2
        center_y = image_size[0] // 2

        pixel_coords = []
        for obs_pos in observed_positions:
            pixel_x = int((obs_pos[0] + center_x * n) / n)
            pixel_y = int((obs_pos[1] + center_y * n) / n)
            pixel_coords.append((pixel_x, pixel_y))

            try:
                fits_magnification = Mag[pixel_y, pixel_x]
            except IndexError:
                fits_magnification = np.nan
            fits_magnifications.append(fits_magnification)

        # Calculate the predicted flux ratios --------- ********** 
        green_image_index = 1
        predicted_flux_ratios = [magnification / predicted_magnifications[green_image_index] for magnification in predicted_magnifications]

        abs_values = [1 / abs(value) if value != 0 else np.nan for value in fits_magnifications]
        predicted_flux_ratios_observed = [abs_values[i] / abs_values[green_image_index] for i in range(len(abs_values))]

        # Print results for the current model
        output_file.write(f"\nModel: {output_filename}\n")
        print(f"\nModel: {output_filename}\n")

        for i, (obs_pos, pred_pos) in enumerate(paired_positions):
            pred_sign = "(pred_negative)" if paired_original_magnifications[i] < 0 else "(pred_positive)"
            obs_sign = "(obs_negative)" if fits_magnifications[i] < 0 else "(obs_positive)"

            pred_phrase = "(over_pred)" if predicted_flux_ratios[i] > observed_flux_ratios[i] else "(under_pred)" if predicted_flux_ratios[i] < observed_flux_ratios[i] else ""
            obs_phrase = "(over_obs)" if predicted_flux_ratios_observed[i] > observed_flux_ratios[i] else "(under_obs)" if predicted_flux_ratios_observed[i] < observed_flux_ratios[i] else ""

            pred_tag = "[YES_pred]" if (pred_sign == "(pred_negative)" and pred_phrase == "(over_pred)") or (pred_sign == "(pred_positive)" and pred_phrase == "(under_pred)") else "[NO_pred]" if (pred_sign == "(pred_negative)" and pred_phrase == "(under_pred)") or (pred_sign == "(pred_positive)" and pred_phrase == "(over_pred)") else ""
            obs_tag = "[YES_obs]" if (obs_sign == "(obs_negative)" and obs_phrase == "(over_obs)") or (obs_sign == "(obs_positive)" and obs_phrase == "(under_obs)") else "[NO_obs]" if (obs_sign == "(obs_negative)" and obs_phrase == "(under_obs)") or (obs_sign == "(obs_positive)" and obs_phrase == "(over_obs)") else ""

            result_str = (
                f"Pair {i+1} {pred_sign} {obs_sign} {pred_phrase} {obs_phrase} {pred_tag} {obs_tag}:\n"
                f"Observed Position: {obs_pos}\n"
                f"Closest Predicted Position: {pred_pos}\n"
                f"Observed Flux Ratio: {observed_flux_ratios[i]:.2f}\n"
                f"Predicted Magnification: {paired_original_magnifications[i]:.2f}\n"
                f"Magnification from FITS file: {fits_magnifications[i]:.2f}\n"
                f"Predicted Flux Ratio at Predicted Position: {predicted_flux_ratios[i]:.2f}\n"
                f"Predicted Flux Ratio at Observed Position (FITS): {predicted_flux_ratios_observed[i]:.2f}\n\n"
            )

            output_file.write(result_str)
            print(result_str)

            # Update counters for statistical analysis
            if pred_tag:
                if pred_sign == "(pred_negative)":
                    negative_parity_pred_total += 1
                    if pred_tag == "[YES_pred]":
                        negative_parity_pred_yes += 1
                elif pred_sign == "(pred_positive)":
                    positive_parity_pred_total += 1
                    if pred_tag == "[YES_pred]":
                        positive_parity_pred_yes += 1

            if obs_tag:
                if obs_sign == "(obs_negative)":
                    negative_parity_obs_total += 1
                    if obs_tag == "[YES_obs]":
                        negative_parity_obs_yes += 1
                elif obs_sign == "(obs_positive)":
                    positive_parity_obs_total += 1
                    if obs_tag == "[YES_obs]":
                        positive_parity_obs_yes += 1

            if pred_tag == "[YES_pred]":
                yes_pred_count += 1
            elif pred_tag == "[NO_pred]":
                no_pred_count += 1

            if obs_tag == "[YES_obs]":
                yes_obs_count += 1
            elif obs_tag == "[NO_obs]":
                no_obs_count += 1

# Calculate statistical analysis results
def calculate_percentage(numerator, denominator):
    return (numerator / denominator * 100) if denominator > 0 else 0

# Predicted positions
negative_parity_pred_benefit = calculate_percentage(
    negative_parity_pred_yes,
    negative_parity_pred_total
)

positive_parity_pred_benefit = calculate_percentage(
    positive_parity_pred_yes,
    positive_parity_pred_total
)

# Observed positions
negative_parity_obs_benefit = calculate_percentage(
    negative_parity_obs_yes,
    negative_parity_obs_total
)

positive_parity_obs_benefit = calculate_percentage(
    positive_parity_obs_yes,
    positive_parity_obs_total
)

# Overall benefit
overall_pred_benefit = calculate_percentage(
    yes_pred_count,
    yes_pred_count + no_pred_count
)

overall_obs_benefit = calculate_percentage(
    yes_obs_count,
    yes_obs_count + no_obs_count
)

# Output the statistical results
with open(Output_text_file, 'a') as output_file:
    results = (
        "\nStatistical Analysis:\n"
        f"% of negative parity that will benefit from adding subhaloes (predicted positions): {negative_parity_pred_benefit:.2f}%\n"
        f"% of positive parity that will benefit from adding subhaloes (predicted positions): {positive_parity_pred_benefit:.2f}%\n"
        f"% of negative parity that will benefit from adding subhaloes (observed positions): {negative_parity_obs_benefit:.2f}%\n"
        f"% of positive parity that will benefit from adding subhaloes (observed positions): {positive_parity_obs_benefit:.2f}%\n"
        f"% of images benefited from adding subhaloes (predicted positions): {overall_pred_benefit:.2f}%\n"
        f"% of images benefited from adding subhaloes (observed positions): {overall_obs_benefit:.2f}%\n"
    )
    
    output_file.write(results)
    print(results)

