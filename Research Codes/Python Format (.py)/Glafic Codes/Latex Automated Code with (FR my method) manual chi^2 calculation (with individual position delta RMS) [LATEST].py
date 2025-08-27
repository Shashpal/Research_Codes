#!/usr/bin/env python
# coding: utf-8

# In[1]:


# NEW LATEX AUTOMATED CODE WITH MANUAL CHI^2 CALCULATION with individual position delta RMS and flux ratio differences
# RXJ0911 
import numpy as np
import pandas as pd
from astropy.io import fits
import os

# Function to calculate distance
def calculate_distance(x1, y1, x2, y2):
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

# Function to read observed positions, flux ratios, and uncertainties
def read_observed_positions(input_file):
    observed_positions = []
    observed_flux_ratios = []
    position_uncertainties = []
    flux_uncertainties = []
    with open(input_file, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line
        for line in lines:
            columns = line.split()
            if len(columns) >= 5:
                x1 = float(columns[0])
                y1 = float(columns[1])
                flux_ratio = float(columns[2])
                pos_err = float(columns[3])
                flux_err = float(columns[4])
                observed_positions.append((x1, y1))
                observed_flux_ratios.append(flux_ratio)
                position_uncertainties.append(pos_err)
                flux_uncertainties.append(flux_err)
    return observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties

# Function to read predicted positions and magnifications
def read_predicted_positions(output_file):
    predicted_positions = []
    predicted_magnifications = []
    with open(output_file, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line
        for line in lines:
            if line.startswith("#"):
                continue
            columns = line.split()
            if len(columns) >= 3:
                x2 = float(columns[0])
                y2 = float(columns[1])
                M2 = abs(float(columns[2]))
                predicted_positions.append((x2, y2))
                predicted_magnifications.append(M2)
    return predicted_positions, predicted_magnifications

# Function to pair predicted and observed positions
def pair_positions(observed_positions, predicted_positions, predicted_magnifications):
    paired_positions = []
    for observed_position in observed_positions:
        min_distance = float('inf')
        closest_index = None
        for i, predicted_position in enumerate(predicted_positions):
            distance = calculate_distance(*observed_position, *predicted_position)
            if distance < min_distance:
                min_distance = distance
                closest_index = i
        paired_positions.append((observed_position, predicted_positions[closest_index], predicted_magnifications[closest_index]))
    return paired_positions

# Function to calculate manual chi² with proper uncertainties
def calculate_manual_chi2(observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties, 
                         predicted_positions, predicted_magnifications, green_index=1):
    # Pair positions
    paired_positions = pair_positions(observed_positions, predicted_positions, predicted_magnifications)
    
    # Order predicted positions and magnifications based on observed order
    predicted_positions_ordered = [pair[1] for pair in paired_positions]
    predicted_magnifications_ordered = [pair[2] for pair in paired_positions]
    
    # Calculate predicted flux ratios
    predicted_flux_ratios = [m / predicted_magnifications_ordered[green_index] for m in predicted_magnifications_ordered]
    
    chi2 = 0.0
    
    # Position contributions to chi²
    for i, (observed, predicted, _) in enumerate(paired_positions):
        dx = observed[0] - predicted[0]
        dy = observed[1] - predicted[1]
        pos_err = position_uncertainties[i]
        chi2 += (dx/pos_err)**2 + (dy/pos_err)**2
    
    # Flux ratio contributions to chi² (excluding the normalization image)
    for i, (pred_flux, obs_flux) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)):
        if i != green_index:  # Exclude the green image used for normalization
            flux_err = flux_uncertainties[i]
            chi2 += ((pred_flux - obs_flux)/flux_err)**2
    
    return chi2

# Function to extract the final chi^2 from optresult.dat file
def extract_glafic_chi2(optresult_path):
    chi2_value = None
    with open(optresult_path, 'r') as file:
        lines = file.readlines()
    
    # Find the last occurrence of chi^2 (final/best result)
    for line in lines:
        if line.startswith("chi^2"):
            chi2_value = float(line.split()[2])
    
    return chi2_value

# Function to process file group
def process_file_group(output_file, fits_file, input_file, n):
    observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties = read_observed_positions(input_file)
    predicted_positions, predicted_magnifications = read_predicted_positions(output_file)
    
    paired_positions = pair_positions(observed_positions, predicted_positions, predicted_magnifications)
    
    # Calculate Delta RMS of Image Positions correctly
    squared_distances = 0.0
    individual_distances = []  # Store individual distances for each image
    for observed, predicted, _ in paired_positions:
        distance = calculate_distance(*observed, *predicted)
        individual_distances.append(distance * 1000)  # Convert to mas
        squared_distances += distance ** 2  # Summing squared distances
    
    mean_squared_distances = squared_distances / len(paired_positions)
    delta_rms = np.sqrt(mean_squared_distances)  # Taking square root of the mean of squared distances
    delta_rms *= 1000  # Convert to mas
    
    # Order predicted positions and magnifications based on observed order
    predicted_positions_ordered = [pair[1] for pair in paired_positions]
    predicted_magnifications_ordered = [pair[2] for pair in paired_positions]
    
    # Determine the index of the "Green Image"
    green_index = 1  # CHOOSE WHICH IMAGE TO NORMALISE BY HERE !!!!!!!!!!!!!!!!!!!!!!!

    # Predicted Flux Ratios
    predicted_flux_ratios = [m / predicted_magnifications_ordered[green_index] for m in predicted_magnifications_ordered]

    # Calculate absolute differences between predicted and observed flux ratios (excluding normalization image)
    flux_ratio_differences = []
    for i, (pred_flux, obs_flux) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)):
        if i != green_index:  # Exclude the green image used for normalization
            flux_ratio_differences.append(abs(pred_flux - obs_flux))

    # Calculate manual chi²
    manual_chi2 = calculate_manual_chi2(observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties,
                                       predicted_positions, predicted_magnifications, green_index)

    with fits.open(fits_file) as hdul:
        mag = hdul[0].data[6]
    
    pixel_coords = [(int((x + mag.shape[1] // 2 * n) / n), int((y + mag.shape[0] // 2 * n) / n)) for x, y in observed_positions]
    values = [mag[y, x] for x, y in pixel_coords]
    
    abs_values = [1 / abs(value) for value in values]
    predicted_flux_ratios_observed = [abs_values[i] / abs_values[green_index] for i in range(len(abs_values))]

    # Calculate RMS excluding the green image
    squared_diffs_pred = [(pred - obs) ** 2 for i, (pred, obs) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)) if i != green_index]
    predicted_flux_ratios_rms = np.sqrt(sum(squared_diffs_pred) / len(squared_diffs_pred))

    squared_diffs_obs = [(pred - obs) ** 2 for i, (pred, obs) in enumerate(zip(predicted_flux_ratios_observed, observed_flux_ratios)) if i != green_index]
    predicted_flux_ratios_observed_rms = np.sqrt(sum(squared_diffs_obs) / len(squared_diffs_obs))

    return {
        'delta_rms': delta_rms,
        'individual_position_deltas': individual_distances,  # Delta RMS for each of the 4 positions
        'flux_ratio_differences': flux_ratio_differences,    # Absolute differences for 3 flux ratios
        'predicted_flux_ratios_rms': predicted_flux_ratios_rms,
        'predicted_flux_ratios_observed_rms': predicted_flux_ratios_observed_rms,
        'manual_chi2': manual_chi2
    }

# Function to extract chi^2 and lens models from optresult.dat file (following Code 1's approach exactly)
def extract_chi2_and_lens_profiles(optresult_path):
    chi2_value = None
    profiles_used = []
    sie_count = 0

    with open(optresult_path, 'r') as file:
        lines = file.readlines()

    found_chi2 = False
    for line in lines:
        if line.startswith("chi^2"):
            chi2_value = float(line.split()[2])
            found_chi2 = True
            profiles_used = []
            sie_count = 0
        elif found_chi2:
            if line.startswith("lens   sie"):
                profiles_used.append("SIE")
                sie_count += 1
            elif line.startswith("lens   pow"):
                profiles_used.append("POW")
            elif line.startswith("lens   anfw"):
                profiles_used.append("NFW")
            elif line.startswith("lens   pert"):
                profiles_used.append("Shear")

    if chi2_value is None:
        raise ValueError("No chi^2 value found in the file.")

    profiles_used = list(set(profiles_used))
    return chi2_value, profiles_used, sie_count

# Function to calculate free parameters based on lens profiles (identical to Code 1)
def calculate_free_parameters(profiles_used, sie_count):
    profile_params = {
        "SIE": 5,
        "POW": 6,
        "NFW": 6,
        "Shear": 5
    }
    total_params = sum(profile_params[profile] for profile in profiles_used)
    if sie_count == 2:
        total_params += 5
    return total_params

# Function to calculate BIC (Bayesian Information Criterion)
def calculate_bic(chi2, n_data_points, n_free_params):
    return chi2 + n_free_params * np.log(n_data_points)

# Function to compute BIC using manual chi² but following Code 1's approach for free parameters
def compute_bic_with_tracking(manual_chi2, optresult_path):
    # Use Code 1's approach: extract profiles from optresult.dat file
    glafic_chi2, profiles_used, sie_count = extract_chi2_and_lens_profiles(optresult_path)
    n_free_params = calculate_free_parameters(profiles_used, sie_count)
    
    # Properly adjust n_data_points based on constraint type
    if "PF" in optresult_path:
        n_data_points = 4 * 2 + 3  # positions + flux ratios
    else:
        n_data_points = 4 * 2  # positions only
    
    # Always use manual chi² for BIC calculation
    bic_value = calculate_bic(manual_chi2, n_data_points, n_free_params)
    
    return bic_value, glafic_chi2, n_data_points, n_free_params

# Main function to process multiple files
def process_files(input_file, output_files, fits_files, optresult_files, n):
    results = []
    latex_rows = []
    bics = []  # Store BIC values separately
    glafic_chi2_values = []  # Store glafic chi² for reference
    data_points_used = []  # Store number of data points used
    free_params_used = []  # Store number of free parameters used

    for output_file, fits_file, optresult_file in zip(output_files, fits_files, optresult_files):
        # Process output and fits files
        result = process_file_group(output_file, fits_file, input_file, n)
        results.append(result)
        
        # Calculate BIC using manual chi² but track both values
        bic, glafic_chi2, n_data_points, n_free_params = compute_bic_with_tracking(result['manual_chi2'], optresult_file)
        bics.append(bic)
        glafic_chi2_values.append(glafic_chi2)
        data_points_used.append(n_data_points)
        free_params_used.append(n_free_params)
        
        # Create LaTeX row with additional information
        model_name, constraint_name = determine_model_and_constraint(output_file)
        # Format individual position deltas and flux ratio differences for LaTeX
        pos_deltas_str = " & ".join([f"{delta:.3f}" for delta in result['individual_position_deltas']])
        flux_diffs_str = " & ".join([f"{diff:.3f}" for diff in result['flux_ratio_differences']])
        
        latex_row = f"{model_name} & {constraint_name} & {pos_deltas_str} & {result['delta_rms']:.1f} & {flux_diffs_str} & {result['predicted_flux_ratios_rms']:.2f} & {bic:.1f} \\\\"
        latex_rows.append(latex_row)
    
    # Combine results into a DataFrame
    df_results = pd.DataFrame(results)
    df_results['bic'] = bics  # Add BIC values to the DataFrame
    df_results['glafic_chi2'] = glafic_chi2_values  # Add glafic chi² for reference
    df_results['n_data_points'] = data_points_used  # Add number of data points used
    df_results['n_free_params'] = free_params_used  # Add number of free parameters used
    
    return df_results, latex_rows

# Determine model and constraint names from the filename
def determine_model_and_constraint(file_name):
    model_map = {'S': 'SIE', 'P': 'POW', 'N': 'NFW'}
    
    base_name = os.path.basename(file_name)
    model_key = base_name[3]
    constraint_key = base_name[4:6] if base_name[4:6] in ['PF'] else base_name[4]
    
    model_name = model_map.get(model_key, '')
    if 'R' in base_name:
        model_name += ' + Shear'
    if 'G' in base_name:
        model_name += ' + G2'
    
    constraint_name = 'QSO Pos + FR' if constraint_key == 'PF' else 'QSO Pos'
    
    return model_name, constraint_name

# Example usage with filenames
Input_file = "RXJ0911/SPFC/obs_point.dat"

Output_files = [
    "RXJ0911/SPC/outSP_point.dat",
    "RXJ0911/SPC/outSPR_point.dat",
    "RXJ0911/SPC/outSPG_point.dat",
    "RXJ0911/SPC/outSPGR_point.dat",
    "RXJ0911/SPFC/outSPF_point.dat",
    "RXJ0911/SPFC/outSPFR_point.dat",
    "RXJ0911/SPFC/outSPFG_point.dat",
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
    "RXJ0911/NPFC/outNPFGR_point.dat",
]

Fits_files = [
    "RXJ0911/SPC/outSP_lens.fits",
    "RXJ0911/SPC/outSPR_lens.fits",
    "RXJ0911/SPC/outSPG_lens.fits",
    "RXJ0911/SPC/outSPGR_lens.fits",
    "RXJ0911/SPFC/outSPF_lens.fits",
    "RXJ0911/SPFC/outSPFR_lens.fits",
    "RXJ0911/SPFC/outSPFG_lens.fits",
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

Optresult_files = [
    "RXJ0911/SPC/outSP_optresult.dat",
    "RXJ0911/SPC/outSPR_optresult.dat",
    "RXJ0911/SPC/outSPG_optresult.dat",
    "RXJ0911/SPC/outSPGR_optresult.dat",
    "RXJ0911/SPFC/outSPF_optresult.dat",
    "RXJ0911/SPFC/outSPFR_optresult.dat",
    "RXJ0911/SPFC/outSPFG_optresult.dat",
    "RXJ0911/SPFC/outSPFGR_optresult.dat",
    "RXJ0911/PPC/outPP_optresult.dat",
    "RXJ0911/PPC/outPPR_optresult.dat",
    "RXJ0911/PPC/outPPG_optresult.dat",
    "RXJ0911/PPC/outPPGR_optresult.dat",
    "RXJ0911/PPFC/outPPF_optresult.dat",
    "RXJ0911/PPFC/outPPFR_optresult.dat",
    "RXJ0911/PPFC/outPPFG_optresult.dat",
    "RXJ0911/PPFC/outPPFGR_optresult.dat",
    "RXJ0911/NPC/outNP_optresult.dat",
    "RXJ0911/NPC/outNPR_optresult.dat",
    "RXJ0911/NPC/outNPG_optresult.dat",
    "RXJ0911/NPC/outNPGR_optresult.dat",
    "RXJ0911/NPFC/outNPF_optresult.dat",
    "RXJ0911/NPFC/outNPFR_optresult.dat",
    "RXJ0911/NPFC/outNPFG_optresult.dat",
    "RXJ0911/NPFC/outNPFGR_optresult.dat"
]

n = 0.0025  # Example scale factor
results_df, latex_rows = process_files(Input_file, Output_files, Fits_files, Optresult_files, n)

# Print formatted results
for i, output_file in enumerate(Output_files):
    constraint_type = "PF" if "PF" in Optresult_files[i] else "P"
    print(f"Model: {output_file}")
    print(f"Constraint Type: {constraint_type}")
    print(f"Number of data points used: {results_df.iloc[i]['n_data_points']}")
    print(f"Number of free parameters: {results_df.iloc[i]['n_free_params']}")
    print(f"Delta RMS of Image Positions: {results_df.iloc[i]['delta_rms']:.3f}")
    
    # Print individual position deltas (4 positions)
    print(f"Individual Position Deltas (mas):")
    for j, delta in enumerate(results_df.iloc[i]['individual_position_deltas']):
        print(f"  Image {j+1}: {delta:.3f}")
    
    print(f"Delta RMS of flux ratios at the observed position: {results_df.iloc[i]['predicted_flux_ratios_observed_rms']:.3f}")
    print(f"Delta RMS of flux ratios at the predicted position: {results_df.iloc[i]['predicted_flux_ratios_rms']:.3f}")
    
    # Print absolute flux ratio differences (3 ratios, excluding normalization)
    print(f"Absolute Flux Ratio Differences:")
    flux_diff_indices = [0, 2, 3] if len(results_df.iloc[i]['flux_ratio_differences']) == 3 else range(len(results_df.iloc[i]['flux_ratio_differences']))
    for j, diff in enumerate(results_df.iloc[i]['flux_ratio_differences']):
        img_idx = flux_diff_indices[j] + 1
        print(f"  Image {img_idx}: {diff:.3f}")
    
    print(f"Manual Chi²: {results_df.iloc[i]['manual_chi2']:.3f}")
    if results_df.iloc[i]['glafic_chi2'] is not None:
        print(f"Glafic Chi² (final, reference): {results_df.iloc[i]['glafic_chi2']:.3f}")
        print(f"Chi² Difference (Manual - Glafic): {results_df.iloc[i]['manual_chi2'] - results_df.iloc[i]['glafic_chi2']:.3f}")
    else:
        print(f"Glafic Chi² (reference): Not found")
    print(f"BIC (using manual Chi²): {results_df.iloc[i]['bic']:.3f}")
    
    # Manual verification of BIC calculation
    manual_bic_check = results_df.iloc[i]['manual_chi2'] + results_df.iloc[i]['n_free_params'] * np.log(results_df.iloc[i]['n_data_points'])
    print(f"Manual BIC verification: {results_df.iloc[i]['manual_chi2']:.3f} + {results_df.iloc[i]['n_free_params']} × ln({results_df.iloc[i]['n_data_points']}) = {manual_bic_check:.3f}")
    print()

# Print LaTeX-formatted rows
print("LaTeX Table Rows:")
print("Format: Model & Constraint & Pos1 & Pos2 & Pos3 & Pos4 & Overall ΔRMS & FluxDiff1 & FluxDiff2 & FluxDiff3 & Flux RMS & BIC")
for row in latex_rows:
    print(row)


# In[1]:


# NEW LATEX AUTOMATED CODE WITH MANUAL CHI^2 CALCULATION with individual position delta RMS and flux ratio differences
# PSJ1606 
import numpy as np
import pandas as pd
from astropy.io import fits
import os

# Function to calculate distance
def calculate_distance(x1, y1, x2, y2):
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

# Function to read observed positions, flux ratios, and uncertainties
def read_observed_positions(input_file):
    observed_positions = []
    observed_flux_ratios = []
    position_uncertainties = []
    flux_uncertainties = []
    with open(input_file, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line
        for line in lines:
            columns = line.split()
            if len(columns) >= 5:
                x1 = float(columns[0])
                y1 = float(columns[1])
                flux_ratio = float(columns[2])
                pos_err = float(columns[3])
                flux_err = float(columns[4])
                observed_positions.append((x1, y1))
                observed_flux_ratios.append(flux_ratio)
                position_uncertainties.append(pos_err)
                flux_uncertainties.append(flux_err)
    return observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties

# Function to read predicted positions and magnifications
def read_predicted_positions(output_file):
    predicted_positions = []
    predicted_magnifications = []
    with open(output_file, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line
        for line in lines:
            if line.startswith("#"):
                continue
            columns = line.split()
            if len(columns) >= 3:
                x2 = float(columns[0])
                y2 = float(columns[1])
                M2 = abs(float(columns[2]))
                predicted_positions.append((x2, y2))
                predicted_magnifications.append(M2)
    return predicted_positions, predicted_magnifications

# Function to pair predicted and observed positions
def pair_positions(observed_positions, predicted_positions, predicted_magnifications):
    paired_positions = []
    for observed_position in observed_positions:
        min_distance = float('inf')
        closest_index = None
        for i, predicted_position in enumerate(predicted_positions):
            distance = calculate_distance(*observed_position, *predicted_position)
            if distance < min_distance:
                min_distance = distance
                closest_index = i
        paired_positions.append((observed_position, predicted_positions[closest_index], predicted_magnifications[closest_index]))
    return paired_positions

# Function to calculate manual chi² with proper uncertainties
def calculate_manual_chi2(observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties, 
                         predicted_positions, predicted_magnifications, green_index=1):
    # Pair positions
    paired_positions = pair_positions(observed_positions, predicted_positions, predicted_magnifications)
    
    # Order predicted positions and magnifications based on observed order
    predicted_positions_ordered = [pair[1] for pair in paired_positions]
    predicted_magnifications_ordered = [pair[2] for pair in paired_positions]
    
    # Calculate predicted flux ratios
    predicted_flux_ratios = [m / predicted_magnifications_ordered[green_index] for m in predicted_magnifications_ordered]
    
    chi2 = 0.0
    
    # Position contributions to chi²
    for i, (observed, predicted, _) in enumerate(paired_positions):
        dx = observed[0] - predicted[0]
        dy = observed[1] - predicted[1]
        pos_err = position_uncertainties[i]
        chi2 += (dx/pos_err)**2 + (dy/pos_err)**2
    
    # Flux ratio contributions to chi² (excluding the normalization image)
    for i, (pred_flux, obs_flux) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)):
        if i != green_index:  # Exclude the green image used for normalization
            flux_err = flux_uncertainties[i]
            chi2 += ((pred_flux - obs_flux)/flux_err)**2
    
    return chi2

# Function to extract the final chi^2 from optresult.dat file
def extract_glafic_chi2(optresult_path):
    chi2_value = None
    with open(optresult_path, 'r') as file:
        lines = file.readlines()
    
    # Find the last occurrence of chi^2 (final/best result)
    for line in lines:
        if line.startswith("chi^2"):
            chi2_value = float(line.split()[2])
    
    return chi2_value

# Function to process file group
def process_file_group(output_file, fits_file, input_file, n):
    observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties = read_observed_positions(input_file)
    predicted_positions, predicted_magnifications = read_predicted_positions(output_file)
    
    paired_positions = pair_positions(observed_positions, predicted_positions, predicted_magnifications)
    
    # Calculate Delta RMS of Image Positions correctly
    squared_distances = 0.0
    individual_distances = []  # Store individual distances for each image
    for observed, predicted, _ in paired_positions:
        distance = calculate_distance(*observed, *predicted)
        individual_distances.append(distance * 1000)  # Convert to mas
        squared_distances += distance ** 2  # Summing squared distances
    
    mean_squared_distances = squared_distances / len(paired_positions)
    delta_rms = np.sqrt(mean_squared_distances)  # Taking square root of the mean of squared distances
    delta_rms *= 1000  # Convert to mas
    
    # Order predicted positions and magnifications based on observed order
    predicted_positions_ordered = [pair[1] for pair in paired_positions]
    predicted_magnifications_ordered = [pair[2] for pair in paired_positions]
    
    # Determine the index of the "Green Image"
    green_index = 0  # CHOOSE WHICH IMAGE TO NORMALISE BY HERE !!!!!!!!!!!!!!!!!!!!!!!

    # Predicted Flux Ratios
    predicted_flux_ratios = [m / predicted_magnifications_ordered[green_index] for m in predicted_magnifications_ordered]

    # Calculate absolute differences between predicted and observed flux ratios (excluding normalization image)
    flux_ratio_differences = []
    for i, (pred_flux, obs_flux) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)):
        if i != green_index:  # Exclude the green image used for normalization
            flux_ratio_differences.append(abs(pred_flux - obs_flux))

    # Calculate manual chi²
    manual_chi2 = calculate_manual_chi2(observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties,
                                       predicted_positions, predicted_magnifications, green_index)

    with fits.open(fits_file) as hdul:
        mag = hdul[0].data[6]
    
    pixel_coords = [(int((x + mag.shape[1] // 2 * n) / n), int((y + mag.shape[0] // 2 * n) / n)) for x, y in observed_positions]
    values = [mag[y, x] for x, y in pixel_coords]
    
    abs_values = [1 / abs(value) for value in values]
    predicted_flux_ratios_observed = [abs_values[i] / abs_values[green_index] for i in range(len(abs_values))]

    # Calculate RMS excluding the green image
    squared_diffs_pred = [(pred - obs) ** 2 for i, (pred, obs) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)) if i != green_index]
    predicted_flux_ratios_rms = np.sqrt(sum(squared_diffs_pred) / len(squared_diffs_pred))

    squared_diffs_obs = [(pred - obs) ** 2 for i, (pred, obs) in enumerate(zip(predicted_flux_ratios_observed, observed_flux_ratios)) if i != green_index]
    predicted_flux_ratios_observed_rms = np.sqrt(sum(squared_diffs_obs) / len(squared_diffs_obs))

    return {
        'delta_rms': delta_rms,
        'individual_position_deltas': individual_distances,  # Delta RMS for each of the 4 positions
        'flux_ratio_differences': flux_ratio_differences,    # Absolute differences for 3 flux ratios
        'predicted_flux_ratios_rms': predicted_flux_ratios_rms,
        'predicted_flux_ratios_observed_rms': predicted_flux_ratios_observed_rms,
        'manual_chi2': manual_chi2
    }

# Function to extract chi^2 and lens models from optresult.dat file (following Code 1's approach exactly)
def extract_chi2_and_lens_profiles(optresult_path):
    chi2_value = None
    profiles_used = []
    sie_count = 0

    with open(optresult_path, 'r') as file:
        lines = file.readlines()

    found_chi2 = False
    for line in lines:
        if line.startswith("chi^2"):
            chi2_value = float(line.split()[2])
            found_chi2 = True
            profiles_used = []
            sie_count = 0
        elif found_chi2:
            if line.startswith("lens   sie"):
                profiles_used.append("SIE")
                sie_count += 1
            elif line.startswith("lens   pow"):
                profiles_used.append("POW")
            elif line.startswith("lens   anfw"):
                profiles_used.append("NFW")
            elif line.startswith("lens   pert"):
                profiles_used.append("Shear")

    if chi2_value is None:
        raise ValueError("No chi^2 value found in the file.")

    profiles_used = list(set(profiles_used))
    return chi2_value, profiles_used, sie_count

# Function to calculate free parameters based on lens profiles (identical to Code 1)
def calculate_free_parameters(profiles_used, sie_count):
    profile_params = {
        "SIE": 5,
        "POW": 6,
        "NFW": 6,
        "Shear": 5
    }
    total_params = sum(profile_params[profile] for profile in profiles_used)
    if sie_count == 2:
        total_params += 5
    return total_params

# Function to calculate BIC (Bayesian Information Criterion)
def calculate_bic(chi2, n_data_points, n_free_params):
    return chi2 + n_free_params * np.log(n_data_points)

# Function to compute BIC using manual chi² but following Code 1's approach for free parameters
def compute_bic_with_tracking(manual_chi2, optresult_path):
    # Use Code 1's approach: extract profiles from optresult.dat file
    glafic_chi2, profiles_used, sie_count = extract_chi2_and_lens_profiles(optresult_path)
    n_free_params = calculate_free_parameters(profiles_used, sie_count)
    
    # Properly adjust n_data_points based on constraint type
    if "PF" in optresult_path:
        n_data_points = 4 * 2 + 3  # positions + flux ratios
    else:
        n_data_points = 4 * 2  # positions only
    
    # Always use manual chi² for BIC calculation
    bic_value = calculate_bic(manual_chi2, n_data_points, n_free_params)
    
    return bic_value, glafic_chi2, n_data_points, n_free_params

# Main function to process multiple files
def process_files(input_file, output_files, fits_files, optresult_files, n):
    results = []
    latex_rows = []
    bics = []  # Store BIC values separately
    glafic_chi2_values = []  # Store glafic chi² for reference
    data_points_used = []  # Store number of data points used
    free_params_used = []  # Store number of free parameters used

    for output_file, fits_file, optresult_file in zip(output_files, fits_files, optresult_files):
        # Process output and fits files
        result = process_file_group(output_file, fits_file, input_file, n)
        results.append(result)
        
        # Calculate BIC using manual chi² but track both values
        bic, glafic_chi2, n_data_points, n_free_params = compute_bic_with_tracking(result['manual_chi2'], optresult_file)
        bics.append(bic)
        glafic_chi2_values.append(glafic_chi2)
        data_points_used.append(n_data_points)
        free_params_used.append(n_free_params)
        
        # Create LaTeX row with additional information
        model_name, constraint_name = determine_model_and_constraint(output_file)
        # Format individual position deltas and flux ratio differences for LaTeX
        pos_deltas_str = " & ".join([f"{delta:.3f}" for delta in result['individual_position_deltas']])
        flux_diffs_str = " & ".join([f"{diff:.3f}" for diff in result['flux_ratio_differences']])
        
        latex_row = f"{model_name} & {constraint_name} & {pos_deltas_str} & {result['delta_rms']:.1f} & {flux_diffs_str} & {result['predicted_flux_ratios_rms']:.2f} & {bic:.1f} \\\\"
        latex_rows.append(latex_row)
    
    # Combine results into a DataFrame
    df_results = pd.DataFrame(results)
    df_results['bic'] = bics  # Add BIC values to the DataFrame
    df_results['glafic_chi2'] = glafic_chi2_values  # Add glafic chi² for reference
    df_results['n_data_points'] = data_points_used  # Add number of data points used
    df_results['n_free_params'] = free_params_used  # Add number of free parameters used
    
    return df_results, latex_rows

# Determine model and constraint names from the filename
def determine_model_and_constraint(file_name):
    model_map = {'S': 'SIE', 'P': 'POW', 'N': 'NFW'}
    
    base_name = os.path.basename(file_name)
    model_key = base_name[3]
    constraint_key = base_name[4:6] if base_name[4:6] in ['PF'] else base_name[4]
    
    model_name = model_map.get(model_key, '')
    if 'R' in base_name:
        model_name += ' + Shear'
    if 'G' in base_name:
        model_name += ' + G2'
    
    constraint_name = 'QSO Pos + FR' if constraint_key == 'PF' else 'QSO Pos'
    
    return model_name, constraint_name

# Example usage with filenames
Input_file = "PSJ1606/SPFC/Lobs_point.dat"

Output_files = [
    "PSJ1606/SPC/outSP_point.dat",
    "PSJ1606/SPC/outSPR_point.dat",
    "PSJ1606/SPC/outSPG_point.dat",
    "PSJ1606/SPC/outSPGR_point.dat",
    "PSJ1606/SPFC/outSPF_point.dat",
    "PSJ1606/SPFC/outSPFR_point.dat",
    "PSJ1606/SPFC/outSPFG_point.dat",
    "PSJ1606/SPFC/outSPFGR_point.dat",
    "PSJ1606/PPC/outPP_point.dat",
    "PSJ1606/PPC/outPPR_point.dat",
    "PSJ1606/PPC/outPPG_point.dat",
    "PSJ1606/PPC/outPPGR_point.dat",
    "PSJ1606/PPFC/outPPF_point.dat",
    "PSJ1606/PPFC/outPPFR_point.dat",
    "PSJ1606/PPFC/outPPFG_point.dat",
    "PSJ1606/PPFC/outPPFGR_point.dat",
    "PSJ1606/NPC/outNP_point.dat",
    "PSJ1606/NPC/outNPR_point.dat",
    "PSJ1606/NPC/outNPG_point.dat",
    "PSJ1606/NPC/outNPGR_point.dat",
    "PSJ1606/NPFC/outNPF_point.dat",
    "PSJ1606/NPFC/outNPFR_point.dat",
    "PSJ1606/NPFC/outNPFG_point.dat",
    "PSJ1606/NPFC/outNPFGR_point.dat"
]

Fits_files = [
    "PSJ1606/SPC/outSP_lens.fits",
    "PSJ1606/SPC/outSPR_lens.fits",
    "PSJ1606/SPC/outSPG_lens.fits",
    "PSJ1606/SPC/outSPGR_lens.fits",
    "PSJ1606/SPFC/outSPF_lens.fits",
    "PSJ1606/SPFC/outSPFR_lens.fits",
    "PSJ1606/SPFC/outSPFG_lens.fits",
    "PSJ1606/SPFC/outSPFGR_lens.fits",
    "PSJ1606/PPC/outPP_lens.fits",
    "PSJ1606/PPC/outPPR_lens.fits",
    "PSJ1606/PPC/outPPG_lens.fits",
    "PSJ1606/PPC/outPPGR_lens.fits",
    "PSJ1606/PPFC/outPPF_lens.fits",
    "PSJ1606/PPFC/outPPFR_lens.fits",
    "PSJ1606/PPFC/outPPFG_lens.fits",
    "PSJ1606/PPFC/outPPFGR_lens.fits",
    "PSJ1606/NPC/outNP_lens.fits",
    "PSJ1606/NPC/outNPR_lens.fits",
    "PSJ1606/NPC/outNPG_lens.fits",
    "PSJ1606/NPC/outNPGR_lens.fits",
    "PSJ1606/NPFC/outNPF_lens.fits",
    "PSJ1606/NPFC/outNPFR_lens.fits",
    "PSJ1606/NPFC/outNPFG_lens.fits",
    "PSJ1606/NPFC/outNPFGR_lens.fits"
]

Optresult_files = [
    "PSJ1606/SPC/outSP_optresult.dat",
    "PSJ1606/SPC/outSPR_optresult.dat",
    "PSJ1606/SPC/outSPG_optresult.dat",
    "PSJ1606/SPC/outSPGR_optresult.dat",
    "PSJ1606/SPFC/outSPF_optresult.dat",
    "PSJ1606/SPFC/outSPFR_optresult.dat",
    "PSJ1606/SPFC/outSPFG_optresult.dat",
    "PSJ1606/SPFC/outSPFGR_optresult.dat",
    "PSJ1606/PPC/outPP_optresult.dat",
    "PSJ1606/PPC/outPPR_optresult.dat",
    "PSJ1606/PPC/outPPG_optresult.dat",
    "PSJ1606/PPC/outPPGR_optresult.dat",
    "PSJ1606/PPFC/outPPF_optresult.dat",
    "PSJ1606/PPFC/outPPFR_optresult.dat",
    "PSJ1606/PPFC/outPPFG_optresult.dat",
    "PSJ1606/PPFC/outPPFGR_optresult.dat",
    "PSJ1606/NPC/outNP_optresult.dat",
    "PSJ1606/NPC/outNPR_optresult.dat",
    "PSJ1606/NPC/outNPG_optresult.dat",
    "PSJ1606/NPC/outNPGR_optresult.dat",
    "PSJ1606/NPFC/outNPF_optresult.dat",
    "PSJ1606/NPFC/outNPFR_optresult.dat",
    "PSJ1606/NPFC/outNPFG_optresult.dat",
    "PSJ1606/NPFC/outNPFGR_optresult.dat"
]

n = 0.0025  # Example scale factor
results_df, latex_rows = process_files(Input_file, Output_files, Fits_files, Optresult_files, n)

# Print formatted results
for i, output_file in enumerate(Output_files):
    constraint_type = "PF" if "PF" in Optresult_files[i] else "P"
    print(f"Model: {output_file}")
    print(f"Constraint Type: {constraint_type}")
    print(f"Number of data points used: {results_df.iloc[i]['n_data_points']}")
    print(f"Number of free parameters: {results_df.iloc[i]['n_free_params']}")
    print(f"Delta RMS of Image Positions: {results_df.iloc[i]['delta_rms']:.3f}")
    
    # Print individual position deltas (4 positions)
    print(f"Individual Position Deltas (mas):")
    for j, delta in enumerate(results_df.iloc[i]['individual_position_deltas']):
        print(f"  Image {j+1}: {delta:.3f}")
    
    print(f"Delta RMS of flux ratios at the observed position: {results_df.iloc[i]['predicted_flux_ratios_observed_rms']:.3f}")
    print(f"Delta RMS of flux ratios at the predicted position: {results_df.iloc[i]['predicted_flux_ratios_rms']:.3f}")
    
    # Print absolute flux ratio differences (3 ratios, excluding normalization)
    print(f"Absolute Flux Ratio Differences:")
    flux_diff_indices = [0, 2, 3] if len(results_df.iloc[i]['flux_ratio_differences']) == 3 else range(len(results_df.iloc[i]['flux_ratio_differences']))
    for j, diff in enumerate(results_df.iloc[i]['flux_ratio_differences']):
        img_idx = flux_diff_indices[j] + 1
        print(f"  Image {img_idx}: {diff:.3f}")
    
    print(f"Manual Chi²: {results_df.iloc[i]['manual_chi2']:.3f}")
    if results_df.iloc[i]['glafic_chi2'] is not None:
        print(f"Glafic Chi² (final, reference): {results_df.iloc[i]['glafic_chi2']:.3f}")
        print(f"Chi² Difference (Manual - Glafic): {results_df.iloc[i]['manual_chi2'] - results_df.iloc[i]['glafic_chi2']:.3f}")
    else:
        print(f"Glafic Chi² (reference): Not found")
    print(f"BIC (using manual Chi²): {results_df.iloc[i]['bic']:.3f}")
    
    # Manual verification of BIC calculation
    manual_bic_check = results_df.iloc[i]['manual_chi2'] + results_df.iloc[i]['n_free_params'] * np.log(results_df.iloc[i]['n_data_points'])
    print(f"Manual BIC verification: {results_df.iloc[i]['manual_chi2']:.3f} + {results_df.iloc[i]['n_free_params']} × ln({results_df.iloc[i]['n_data_points']}) = {manual_bic_check:.3f}")
    print()

# Print LaTeX-formatted rows
print("LaTeX Table Rows:")
print("Format: Model & Constraint & Pos1 & Pos2 & Pos3 & Pos4 & Overall ΔRMS & FluxDiff1 & FluxDiff2 & FluxDiff3 & Flux RMS & BIC")
for row in latex_rows:
    print(row)


# In[3]:


# NEW LATEX AUTOMATED CODE WITH MANUAL CHI^2 CALCULATION with individual position delta RMS and flux ratio differences
# WFI2033
import numpy as np
import pandas as pd
from astropy.io import fits
import os

# Function to calculate distance
def calculate_distance(x1, y1, x2, y2):
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

# Function to read observed positions, flux ratios, and uncertainties
def read_observed_positions(input_file):
    observed_positions = []
    observed_flux_ratios = []
    position_uncertainties = []
    flux_uncertainties = []
    with open(input_file, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line
        for line in lines:
            columns = line.split()
            if len(columns) >= 5:
                x1 = float(columns[0])
                y1 = float(columns[1])
                flux_ratio = float(columns[2])
                pos_err = float(columns[3])
                flux_err = float(columns[4])
                observed_positions.append((x1, y1))
                observed_flux_ratios.append(flux_ratio)
                position_uncertainties.append(pos_err)
                flux_uncertainties.append(flux_err)
    return observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties

# Function to read predicted positions and magnifications
def read_predicted_positions(output_file):
    predicted_positions = []
    predicted_magnifications = []
    with open(output_file, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line
        for line in lines:
            if line.startswith("#"):
                continue
            columns = line.split()
            if len(columns) >= 3:
                x2 = float(columns[0])
                y2 = float(columns[1])
                M2 = abs(float(columns[2]))
                predicted_positions.append((x2, y2))
                predicted_magnifications.append(M2)
    return predicted_positions, predicted_magnifications

# Function to pair predicted and observed positions
def pair_positions(observed_positions, predicted_positions, predicted_magnifications):
    paired_positions = []
    for observed_position in observed_positions:
        min_distance = float('inf')
        closest_index = None
        for i, predicted_position in enumerate(predicted_positions):
            distance = calculate_distance(*observed_position, *predicted_position)
            if distance < min_distance:
                min_distance = distance
                closest_index = i
        paired_positions.append((observed_position, predicted_positions[closest_index], predicted_magnifications[closest_index]))
    return paired_positions

# Function to calculate manual chi² with proper uncertainties
def calculate_manual_chi2(observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties, 
                         predicted_positions, predicted_magnifications, green_index=1):
    # Pair positions
    paired_positions = pair_positions(observed_positions, predicted_positions, predicted_magnifications)
    
    # Order predicted positions and magnifications based on observed order
    predicted_positions_ordered = [pair[1] for pair in paired_positions]
    predicted_magnifications_ordered = [pair[2] for pair in paired_positions]
    
    # Calculate predicted flux ratios
    predicted_flux_ratios = [m / predicted_magnifications_ordered[green_index] for m in predicted_magnifications_ordered]
    
    chi2 = 0.0
    
    # Position contributions to chi²
    for i, (observed, predicted, _) in enumerate(paired_positions):
        dx = observed[0] - predicted[0]
        dy = observed[1] - predicted[1]
        pos_err = position_uncertainties[i]
        chi2 += (dx/pos_err)**2 + (dy/pos_err)**2
    
    # Flux ratio contributions to chi² (excluding the normalization image)
    for i, (pred_flux, obs_flux) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)):
        if i != green_index:  # Exclude the green image used for normalization
            flux_err = flux_uncertainties[i]
            chi2 += ((pred_flux - obs_flux)/flux_err)**2
    
    return chi2

# Function to extract the final chi^2 from optresult.dat file
def extract_glafic_chi2(optresult_path):
    chi2_value = None
    with open(optresult_path, 'r') as file:
        lines = file.readlines()
    
    # Find the last occurrence of chi^2 (final/best result)
    for line in lines:
        if line.startswith("chi^2"):
            chi2_value = float(line.split()[2])
    
    return chi2_value

# Function to process file group
def process_file_group(output_file, fits_file, input_file, n):
    observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties = read_observed_positions(input_file)
    predicted_positions, predicted_magnifications = read_predicted_positions(output_file)
    
    paired_positions = pair_positions(observed_positions, predicted_positions, predicted_magnifications)
    
    # Calculate Delta RMS of Image Positions correctly
    squared_distances = 0.0
    individual_distances = []  # Store individual distances for each image
    for observed, predicted, _ in paired_positions:
        distance = calculate_distance(*observed, *predicted)
        individual_distances.append(distance * 1000)  # Convert to mas
        squared_distances += distance ** 2  # Summing squared distances
    
    mean_squared_distances = squared_distances / len(paired_positions)
    delta_rms = np.sqrt(mean_squared_distances)  # Taking square root of the mean of squared distances
    delta_rms *= 1000  # Convert to mas
    
    # Order predicted positions and magnifications based on observed order
    predicted_positions_ordered = [pair[1] for pair in paired_positions]
    predicted_magnifications_ordered = [pair[2] for pair in paired_positions]
    
    # Determine the index of the "Green Image"
    green_index = 0  # CHOOSE WHICH IMAGE TO NORMALISE BY HERE !!!!!!!!!!!!!!!!!!!!!!!

    # Predicted Flux Ratios
    predicted_flux_ratios = [m / predicted_magnifications_ordered[green_index] for m in predicted_magnifications_ordered]

    # Calculate absolute differences between predicted and observed flux ratios (excluding normalization image)
    flux_ratio_differences = []
    for i, (pred_flux, obs_flux) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)):
        if i != green_index:  # Exclude the green image used for normalization
            flux_ratio_differences.append(abs(pred_flux - obs_flux))

    # Calculate manual chi²
    manual_chi2 = calculate_manual_chi2(observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties,
                                       predicted_positions, predicted_magnifications, green_index)

    with fits.open(fits_file) as hdul:
        mag = hdul[0].data[6]
    
    pixel_coords = [(int((x + mag.shape[1] // 2 * n) / n), int((y + mag.shape[0] // 2 * n) / n)) for x, y in observed_positions]
    values = [mag[y, x] for x, y in pixel_coords]
    
    abs_values = [1 / abs(value) for value in values]
    predicted_flux_ratios_observed = [abs_values[i] / abs_values[green_index] for i in range(len(abs_values))]

    # Calculate RMS excluding the green image
    squared_diffs_pred = [(pred - obs) ** 2 for i, (pred, obs) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)) if i != green_index]
    predicted_flux_ratios_rms = np.sqrt(sum(squared_diffs_pred) / len(squared_diffs_pred))

    squared_diffs_obs = [(pred - obs) ** 2 for i, (pred, obs) in enumerate(zip(predicted_flux_ratios_observed, observed_flux_ratios)) if i != green_index]
    predicted_flux_ratios_observed_rms = np.sqrt(sum(squared_diffs_obs) / len(squared_diffs_obs))

    return {
        'delta_rms': delta_rms,
        'individual_position_deltas': individual_distances,  # Delta RMS for each of the 4 positions
        'flux_ratio_differences': flux_ratio_differences,    # Absolute differences for 3 flux ratios
        'predicted_flux_ratios_rms': predicted_flux_ratios_rms,
        'predicted_flux_ratios_observed_rms': predicted_flux_ratios_observed_rms,
        'manual_chi2': manual_chi2
    }

# Function to extract chi^2 and lens models from optresult.dat file (following Code 1's approach exactly)
def extract_chi2_and_lens_profiles(optresult_path):
    chi2_value = None
    profiles_used = []
    sie_count = 0

    with open(optresult_path, 'r') as file:
        lines = file.readlines()

    found_chi2 = False
    for line in lines:
        if line.startswith("chi^2"):
            chi2_value = float(line.split()[2])
            found_chi2 = True
            profiles_used = []
            sie_count = 0
        elif found_chi2:
            if line.startswith("lens   sie"):
                profiles_used.append("SIE")
                sie_count += 1
            elif line.startswith("lens   pow"):
                profiles_used.append("POW")
            elif line.startswith("lens   anfw"):
                profiles_used.append("NFW")
            elif line.startswith("lens   pert"):
                profiles_used.append("Shear")

    if chi2_value is None:
        raise ValueError("No chi^2 value found in the file.")

    profiles_used = list(set(profiles_used))
    return chi2_value, profiles_used, sie_count

# Function to calculate free parameters based on lens profiles (identical to Code 1)
def calculate_free_parameters(profiles_used, sie_count):
    profile_params = {
        "SIE": 5,
        "POW": 6,
        "NFW": 6,
        "Shear": 5
    }
    total_params = sum(profile_params[profile] for profile in profiles_used)
    if sie_count == 2:
        total_params += 5
    return total_params

# Function to calculate BIC (Bayesian Information Criterion)
def calculate_bic(chi2, n_data_points, n_free_params):
    return chi2 + n_free_params * np.log(n_data_points)

# Function to compute BIC using manual chi² but following Code 1's approach for free parameters
def compute_bic_with_tracking(manual_chi2, optresult_path):
    # Use Code 1's approach: extract profiles from optresult.dat file
    glafic_chi2, profiles_used, sie_count = extract_chi2_and_lens_profiles(optresult_path)
    n_free_params = calculate_free_parameters(profiles_used, sie_count)
    
    # Properly adjust n_data_points based on constraint type
    if "PF" in optresult_path:
        n_data_points = 4 * 2 + 3  # positions + flux ratios
    else:
        n_data_points = 4 * 2  # positions only
    
    # Always use manual chi² for BIC calculation
    bic_value = calculate_bic(manual_chi2, n_data_points, n_free_params)
    
    return bic_value, glafic_chi2, n_data_points, n_free_params

# Main function to process multiple files
def process_files(input_file, output_files, fits_files, optresult_files, n):
    results = []
    latex_rows = []
    bics = []  # Store BIC values separately
    glafic_chi2_values = []  # Store glafic chi² for reference
    data_points_used = []  # Store number of data points used
    free_params_used = []  # Store number of free parameters used

    for output_file, fits_file, optresult_file in zip(output_files, fits_files, optresult_files):
        # Process output and fits files
        result = process_file_group(output_file, fits_file, input_file, n)
        results.append(result)
        
        # Calculate BIC using manual chi² but track both values
        bic, glafic_chi2, n_data_points, n_free_params = compute_bic_with_tracking(result['manual_chi2'], optresult_file)
        bics.append(bic)
        glafic_chi2_values.append(glafic_chi2)
        data_points_used.append(n_data_points)
        free_params_used.append(n_free_params)
        
        # Create LaTeX row with additional information
        model_name, constraint_name = determine_model_and_constraint(output_file)
        # Format individual position deltas and flux ratio differences for LaTeX
        pos_deltas_str = " & ".join([f"{delta:.3f}" for delta in result['individual_position_deltas']])
        flux_diffs_str = " & ".join([f"{diff:.3f}" for diff in result['flux_ratio_differences']])
        
        latex_row = f"{model_name} & {constraint_name} & {pos_deltas_str} & {result['delta_rms']:.1f} & {flux_diffs_str} & {result['predicted_flux_ratios_rms']:.2f} & {bic:.1f} \\\\"
        latex_rows.append(latex_row)
    
    # Combine results into a DataFrame
    df_results = pd.DataFrame(results)
    df_results['bic'] = bics  # Add BIC values to the DataFrame
    df_results['glafic_chi2'] = glafic_chi2_values  # Add glafic chi² for reference
    df_results['n_data_points'] = data_points_used  # Add number of data points used
    df_results['n_free_params'] = free_params_used  # Add number of free parameters used
    
    return df_results, latex_rows

# Determine model and constraint names from the filename
def determine_model_and_constraint(file_name):
    model_map = {'S': 'SIE', 'P': 'POW', 'N': 'NFW'}
    
    base_name = os.path.basename(file_name)
    model_key = base_name[3]
    constraint_key = base_name[4:6] if base_name[4:6] in ['PF'] else base_name[4]
    
    model_name = model_map.get(model_key, '')
    if 'R' in base_name:
        model_name += ' + Shear'
    if 'G' in base_name:
        model_name += ' + G2'
    
    constraint_name = 'QSO Pos + FR' if constraint_key == 'PF' else 'QSO Pos'
    
    return model_name, constraint_name

# Example usage with filenames
Input_file = "WFI2033/obs_point.dat"

Output_files = [
    "WFI2033/SP/outSP_point.dat",
    "WFI2033/SPR/outSPR_point.dat",
    "WFI2033/SPG/outSPG_point.dat",
    "WFI2033/SPGR/outSPGR_point.dat",
    "WFI2033/SPF/outSPF_point.dat",
    "WFI2033/SPFR/outSPFR_point.dat",
    "WFI2033/SPFG/outSPFG_point.dat",
    "WFI2033/SPFGR/outSPFGR_point.dat",
    "WFI2033/PP/outPP_point.dat",
    "WFI2033/PPR/outPPR_point.dat",
    "WFI2033/PPG/outPPG_point.dat",
    "WFI2033/PPGR/outPPGR_point.dat",
    "WFI2033/PPF/outPPF_point.dat",
    "WFI2033/PPFR/outPPFR_point.dat",
    "WFI2033/PPFG/outPPFG_point.dat",
    "WFI2033/PPFGR/outPPFGR_point.dat",
    "WFI2033/NP/outNP_point.dat",
    "WFI2033/NPR/outNPR_point.dat",
    "WFI2033/NPG/outNPG_point.dat",
    "WFI2033/NPGR/outNPGR_point.dat",
    "WFI2033/NPF/outNPF_point.dat",
    "WFI2033/NPFR/outNPFR_point.dat",
    "WFI2033/NPFG/outNPFG_point.dat",
    "WFI2033/NPFGR/outNPFGR_point.dat"
]

Fits_files = [
    "WFI2033/SP/outSP_lens.fits",
    "WFI2033/SPR/outSPR_lens.fits",
    "WFI2033/SPG/outSPG_lens.fits",
    "WFI2033/SPGR/outSPGR_lens.fits",
    "WFI2033/SPF/outSPF_lens.fits",
    "WFI2033/SPFR/outSPFR_lens.fits",
    "WFI2033/SPFG/outSPFG_lens.fits",
    "WFI2033/SPFGR/outSPFGR_lens.fits",
    "WFI2033/PP/outPP_lens.fits",
    "WFI2033/PPR/outPPR_lens.fits",
    "WFI2033/PPG/outPPG_lens.fits",
    "WFI2033/PPGR/outPPGR_lens.fits",
    "WFI2033/PPF/outPPF_lens.fits",
    "WFI2033/PPFR/outPPFR_lens.fits",
    "WFI2033/PPFG/outPPFG_lens.fits",
    "WFI2033/PPFGR/outPPFGR_lens.fits",
    "WFI2033/NP/outNP_lens.fits",
    "WFI2033/NPR/outNPR_lens.fits",
    "WFI2033/NPG/outNPG_lens.fits",
    "WFI2033/NPGR/outNPGR_lens.fits",
    "WFI2033/NPF/outNPF_lens.fits",
    "WFI2033/NPFR/outNPFR_lens.fits",
    "WFI2033/NPFG/outNPFG_lens.fits",
    "WFI2033/NPFGR/outNPFGR_lens.fits"
]

Optresult_files = [
    "WFI2033/SP/outSP_optresult.dat",
    "WFI2033/SPR/outSPR_optresult.dat",
    "WFI2033/SPG/outSPG_optresult.dat",
    "WFI2033/SPGR/outSPGR_optresult.dat",
    "WFI2033/SPF/outSPF_optresult.dat",
    "WFI2033/SPFR/outSPFR_optresult.dat",
    "WFI2033/SPFG/outSPFG_optresult.dat",
    "WFI2033/SPFGR/outSPFGR_optresult.dat",
    "WFI2033/PP/outPP_optresult.dat",
    "WFI2033/PPR/outPPR_optresult.dat",
    "WFI2033/PPG/outPPG_optresult.dat",
    "WFI2033/PPGR/outPPGR_optresult.dat",
    "WFI2033/PPF/outPPF_optresult.dat",
    "WFI2033/PPFR/outPPFR_optresult.dat",
    "WFI2033/PPFG/outPPFG_optresult.dat",
    "WFI2033/PPFGR/outPPFGR_optresult.dat",
    "WFI2033/NP/outNP_optresult.dat",
    "WFI2033/NPR/outNPR_optresult.dat",
    "WFI2033/NPG/outNPG_optresult.dat",
    "WFI2033/NPGR/outNPGR_optresult.dat",
    "WFI2033/NPF/outNPF_optresult.dat",
    "WFI2033/NPFR/outNPFR_optresult.dat",
    "WFI2033/NPFG/outNPFG_optresult.dat",
    "WFI2033/NPFGR/outNPFGR_optresult.dat"
]

n = 0.0025  # Example scale factor
results_df, latex_rows = process_files(Input_file, Output_files, Fits_files, Optresult_files, n)

# Print formatted results
for i, output_file in enumerate(Output_files):
    constraint_type = "PF" if "PF" in Optresult_files[i] else "P"
    print(f"Model: {output_file}")
    print(f"Constraint Type: {constraint_type}")
    print(f"Number of data points used: {results_df.iloc[i]['n_data_points']}")
    print(f"Number of free parameters: {results_df.iloc[i]['n_free_params']}")
    print(f"Delta RMS of Image Positions: {results_df.iloc[i]['delta_rms']:.3f}")
    
    # Print individual position deltas (4 positions)
    print(f"Individual Position Deltas (mas):")
    for j, delta in enumerate(results_df.iloc[i]['individual_position_deltas']):
        print(f"  Image {j+1}: {delta:.3f}")
    
    print(f"Delta RMS of flux ratios at the observed position: {results_df.iloc[i]['predicted_flux_ratios_observed_rms']:.3f}")
    print(f"Delta RMS of flux ratios at the predicted position: {results_df.iloc[i]['predicted_flux_ratios_rms']:.3f}")
    
    # Print absolute flux ratio differences (3 ratios, excluding normalization)
    print(f"Absolute Flux Ratio Differences:")
    flux_diff_indices = [0, 2, 3] if len(results_df.iloc[i]['flux_ratio_differences']) == 3 else range(len(results_df.iloc[i]['flux_ratio_differences']))
    for j, diff in enumerate(results_df.iloc[i]['flux_ratio_differences']):
        img_idx = flux_diff_indices[j] + 1
        print(f"  Image {img_idx}: {diff:.3f}")
    
    print(f"Manual Chi²: {results_df.iloc[i]['manual_chi2']:.3f}")
    if results_df.iloc[i]['glafic_chi2'] is not None:
        print(f"Glafic Chi² (final, reference): {results_df.iloc[i]['glafic_chi2']:.3f}")
        print(f"Chi² Difference (Manual - Glafic): {results_df.iloc[i]['manual_chi2'] - results_df.iloc[i]['glafic_chi2']:.3f}")
    else:
        print(f"Glafic Chi² (reference): Not found")
    print(f"BIC (using manual Chi²): {results_df.iloc[i]['bic']:.3f}")
    
    # Manual verification of BIC calculation
    manual_bic_check = results_df.iloc[i]['manual_chi2'] + results_df.iloc[i]['n_free_params'] * np.log(results_df.iloc[i]['n_data_points'])
    print(f"Manual BIC verification: {results_df.iloc[i]['manual_chi2']:.3f} + {results_df.iloc[i]['n_free_params']} × ln({results_df.iloc[i]['n_data_points']}) = {manual_bic_check:.3f}")
    print()

# Print LaTeX-formatted rows
print("LaTeX Table Rows:")
print("Format: Model & Constraint & Pos1 & Pos2 & Pos3 & Pos4 & Overall ΔRMS & FluxDiff1 & FluxDiff2 & FluxDiff3 & Flux RMS & BIC")
for row in latex_rows:
    print(row)


# In[6]:


# NEW LATEX AUTOMATED CODE WITH MANUAL CHI^2 CALCULATION with individual position delta RMS and flux ratio differences
# WFI2033
import numpy as np
import pandas as pd
from astropy.io import fits
import os

# Function to calculate distance
def calculate_distance(x1, y1, x2, y2):
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

# Function to read observed positions, flux ratios, and uncertainties
def read_observed_positions(input_file):
    observed_positions = []
    observed_flux_ratios = []
    position_uncertainties = []
    flux_uncertainties = []
    with open(input_file, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line
        for line in lines:
            columns = line.split()
            if len(columns) >= 5:
                x1 = float(columns[0])
                y1 = float(columns[1])
                flux_ratio = float(columns[2])
                pos_err = float(columns[3])
                flux_err = float(columns[4])
                observed_positions.append((x1, y1))
                observed_flux_ratios.append(flux_ratio)
                position_uncertainties.append(pos_err)
                flux_uncertainties.append(flux_err)
    return observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties

# Function to read predicted positions and magnifications
def read_predicted_positions(output_file):
    predicted_positions = []
    predicted_magnifications = []
    with open(output_file, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line
        for line in lines:
            if line.startswith("#"):
                continue
            columns = line.split()
            if len(columns) >= 3:
                x2 = float(columns[0])
                y2 = float(columns[1])
                M2 = abs(float(columns[2]))
                predicted_positions.append((x2, y2))
                predicted_magnifications.append(M2)
    return predicted_positions, predicted_magnifications

# Function to pair predicted and observed positions
def pair_positions(observed_positions, predicted_positions, predicted_magnifications):
    paired_positions = []
    for observed_position in observed_positions:
        min_distance = float('inf')
        closest_index = None
        for i, predicted_position in enumerate(predicted_positions):
            distance = calculate_distance(*observed_position, *predicted_position)
            if distance < min_distance:
                min_distance = distance
                closest_index = i
        paired_positions.append((observed_position, predicted_positions[closest_index], predicted_magnifications[closest_index]))
    return paired_positions

# Function to calculate manual chi² with proper uncertainties
def calculate_manual_chi2(observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties, 
                         predicted_positions, predicted_magnifications, green_index=1):
    # Pair positions
    paired_positions = pair_positions(observed_positions, predicted_positions, predicted_magnifications)
    
    # Order predicted positions and magnifications based on observed order
    predicted_positions_ordered = [pair[1] for pair in paired_positions]
    predicted_magnifications_ordered = [pair[2] for pair in paired_positions]
    
    # Calculate predicted flux ratios
    predicted_flux_ratios = [m / predicted_magnifications_ordered[green_index] for m in predicted_magnifications_ordered]
    
    chi2 = 0.0
    
    # Position contributions to chi²
    for i, (observed, predicted, _) in enumerate(paired_positions):
        dx = observed[0] - predicted[0]
        dy = observed[1] - predicted[1]
        pos_err = position_uncertainties[i]
        chi2 += (dx/pos_err)**2 + (dy/pos_err)**2
    
    # Flux ratio contributions to chi² (excluding the normalization image)
    for i, (pred_flux, obs_flux) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)):
        if i != green_index:  # Exclude the green image used for normalization
            flux_err = flux_uncertainties[i]
            chi2 += ((pred_flux - obs_flux)/flux_err)**2
    
    return chi2

# Function to extract the final chi^2 from optresult.dat file
def extract_glafic_chi2(optresult_path):
    chi2_value = None
    with open(optresult_path, 'r') as file:
        lines = file.readlines()
    
    # Find the last occurrence of chi^2 (final/best result)
    for line in lines:
        if line.startswith("chi^2"):
            chi2_value = float(line.split()[2])
    
    return chi2_value

# Function to process file group
def process_file_group(output_file, fits_file, input_file, n):
    observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties = read_observed_positions(input_file)
    predicted_positions, predicted_magnifications = read_predicted_positions(output_file)
    
    paired_positions = pair_positions(observed_positions, predicted_positions, predicted_magnifications)
    
    # Calculate Delta RMS of Image Positions correctly
    squared_distances = 0.0
    individual_distances = []  # Store individual distances for each image
    for observed, predicted, _ in paired_positions:
        distance = calculate_distance(*observed, *predicted)
        individual_distances.append(distance * 1000)  # Convert to mas
        squared_distances += distance ** 2  # Summing squared distances
    
    mean_squared_distances = squared_distances / len(paired_positions)
    delta_rms = np.sqrt(mean_squared_distances)  # Taking square root of the mean of squared distances
    delta_rms *= 1000  # Convert to mas
    
    # Order predicted positions and magnifications based on observed order
    predicted_positions_ordered = [pair[1] for pair in paired_positions]
    predicted_magnifications_ordered = [pair[2] for pair in paired_positions]
    
    # Determine the index of the "Green Image"
    green_index = 0  # CHOOSE WHICH IMAGE TO NORMALISE BY HERE !!!!!!!!!!!!!!!!!!!!!!!

    # Predicted Flux Ratios
    predicted_flux_ratios = [m / predicted_magnifications_ordered[green_index] for m in predicted_magnifications_ordered]

    # Calculate absolute differences between predicted and observed flux ratios (excluding normalization image)
    flux_ratio_differences = []
    for i, (pred_flux, obs_flux) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)):
        if i != green_index:  # Exclude the green image used for normalization
            flux_ratio_differences.append(abs(pred_flux - obs_flux))

    # Calculate manual chi²
    manual_chi2 = calculate_manual_chi2(observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties,
                                       predicted_positions, predicted_magnifications, green_index)

    with fits.open(fits_file) as hdul:
        mag = hdul[0].data[6]
    
    pixel_coords = [(int((x + mag.shape[1] // 2 * n) / n), int((y + mag.shape[0] // 2 * n) / n)) for x, y in observed_positions]
    values = [mag[y, x] for x, y in pixel_coords]
    
    abs_values = [1 / abs(value) for value in values]
    predicted_flux_ratios_observed = [abs_values[i] / abs_values[green_index] for i in range(len(abs_values))]

    # Calculate RMS excluding the green image
    squared_diffs_pred = [(pred - obs) ** 2 for i, (pred, obs) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)) if i != green_index]
    predicted_flux_ratios_rms = np.sqrt(sum(squared_diffs_pred) / len(squared_diffs_pred))

    squared_diffs_obs = [(pred - obs) ** 2 for i, (pred, obs) in enumerate(zip(predicted_flux_ratios_observed, observed_flux_ratios)) if i != green_index]
    predicted_flux_ratios_observed_rms = np.sqrt(sum(squared_diffs_obs) / len(squared_diffs_obs))

    return {
        'delta_rms': delta_rms,
        'individual_position_deltas': individual_distances,  # Delta RMS for each of the 4 positions
        'flux_ratio_differences': flux_ratio_differences,    # Absolute differences for 3 flux ratios
        'predicted_flux_ratios_rms': predicted_flux_ratios_rms,
        'predicted_flux_ratios_observed_rms': predicted_flux_ratios_observed_rms,
        'manual_chi2': manual_chi2
    }

# Function to extract chi^2 and lens models from optresult.dat file (following Code 1's approach exactly)
def extract_chi2_and_lens_profiles(optresult_path):
    chi2_value = None
    profiles_used = []
    sie_count = 0

    with open(optresult_path, 'r') as file:
        lines = file.readlines()

    found_chi2 = False
    for line in lines:
        if line.startswith("chi^2"):
            chi2_value = float(line.split()[2])
            found_chi2 = True
            profiles_used = []
            sie_count = 0
        elif found_chi2:
            if line.startswith("lens   sie"):
                profiles_used.append("SIE")
                sie_count += 1
            elif line.startswith("lens   pow"):
                profiles_used.append("POW")
            elif line.startswith("lens   anfw"):
                profiles_used.append("NFW")
            elif line.startswith("lens   pert"):
                profiles_used.append("Shear")

    if chi2_value is None:
        raise ValueError("No chi^2 value found in the file.")

    profiles_used = list(set(profiles_used))
    return chi2_value, profiles_used, sie_count

# Function to calculate free parameters based on lens profiles (identical to Code 1)
def calculate_free_parameters(profiles_used, sie_count):
    profile_params = {
        "SIE": 5,
        "POW": 6,
        "NFW": 6,
        "Shear": 5
    }
    total_params = sum(profile_params[profile] for profile in profiles_used)
    if sie_count == 2:
        total_params += 5
    return total_params

# Function to calculate BIC (Bayesian Information Criterion)
def calculate_bic(chi2, n_data_points, n_free_params):
    return chi2 + n_free_params * np.log(n_data_points)

# Function to compute BIC using manual chi² but following Code 1's approach for free parameters
def compute_bic_with_tracking(manual_chi2, optresult_path):
    # Use Code 1's approach: extract profiles from optresult.dat file
    glafic_chi2, profiles_used, sie_count = extract_chi2_and_lens_profiles(optresult_path)
    n_free_params = calculate_free_parameters(profiles_used, sie_count)
    
    # Properly adjust n_data_points based on constraint type
    if "PF" in optresult_path:
        n_data_points = 4 * 2 + 3  # positions + flux ratios
    else:
        n_data_points = 4 * 2  # positions only
    
    # Always use manual chi² for BIC calculation
    bic_value = calculate_bic(manual_chi2, n_data_points, n_free_params)
    
    return bic_value, glafic_chi2, n_data_points, n_free_params

# Main function to process multiple files
def process_files(input_file, output_files, fits_files, optresult_files, n):
    results = []
    latex_rows = []
    bics = []  # Store BIC values separately
    glafic_chi2_values = []  # Store glafic chi² for reference
    data_points_used = []  # Store number of data points used
    free_params_used = []  # Store number of free parameters used

    for output_file, fits_file, optresult_file in zip(output_files, fits_files, optresult_files):
        # Process output and fits files
        result = process_file_group(output_file, fits_file, input_file, n)
        results.append(result)
        
        # Calculate BIC using manual chi² but track both values
        bic, glafic_chi2, n_data_points, n_free_params = compute_bic_with_tracking(result['manual_chi2'], optresult_file)
        bics.append(bic)
        glafic_chi2_values.append(glafic_chi2)
        data_points_used.append(n_data_points)
        free_params_used.append(n_free_params)
        
        # Create LaTeX row with additional information
        model_name, constraint_name = determine_model_and_constraint(output_file)
        # Format individual position deltas and flux ratio differences for LaTeX
        pos_deltas_str = " & ".join([f"{delta:.3f}" for delta in result['individual_position_deltas']])
        flux_diffs_str = " & ".join([f"{diff:.3f}" for diff in result['flux_ratio_differences']])
        
        latex_row = f"{model_name} & {constraint_name} & {pos_deltas_str} & {result['delta_rms']:.1f} & {flux_diffs_str} & {result['predicted_flux_ratios_rms']:.2f} & {bic:.1f} \\\\"
        latex_rows.append(latex_row)
    
    # Combine results into a DataFrame
    df_results = pd.DataFrame(results)
    df_results['bic'] = bics  # Add BIC values to the DataFrame
    df_results['glafic_chi2'] = glafic_chi2_values  # Add glafic chi² for reference
    df_results['n_data_points'] = data_points_used  # Add number of data points used
    df_results['n_free_params'] = free_params_used  # Add number of free parameters used
    
    return df_results, latex_rows

# Determine model and constraint names from the filename
def determine_model_and_constraint(file_name):
    model_map = {'S': 'SIE', 'P': 'POW', 'N': 'NFW'}
    
    base_name = os.path.basename(file_name)
    model_key = base_name[3]
    constraint_key = base_name[4:6] if base_name[4:6] in ['PF'] else base_name[4]
    
    model_name = model_map.get(model_key, '')
    if 'R' in base_name:
        model_name += ' + Shear'
    if 'G' in base_name:
        model_name += ' + G2'
    
    constraint_name = 'QSO Pos + FR' if constraint_key == 'PF' else 'QSO Pos'
    
    return model_name, constraint_name

# Example usage with filenames
Input_file = "WFI2033/obs_point.dat"

Output_files = [
    "WFI2033/SP/outSP_point.dat",
    "WFI2033/SPR/outSPR_point.dat",
    "WFI2033/SPG/outSPG_point.dat",
    "WFI2033/SPGR/outSPGR_point.dat",
    "WFI2033/SPF/outSPF_point.dat",
    "WFI2033/SPFR/outSPFR_point.dat",
    "WFI2033/SPFG/outSPFG_point.dat",
    "WFI2033/SPFGR/outSPFGR_point.dat",
    "WFI2033/PP/outPP_point.dat",
    "WFI2033/PPR/outPPR_point.dat",
    "WFI2033/PPG/outPPG_point.dat",
    "WFI2033/PPGR/outPPGR_point.dat",
    "WFI2033/PPF/outPPF_point.dat",
    "WFI2033/PPFR/outPPFR_point.dat",
    "WFI2033/PPFG/outPPFG_point.dat",
    "WFI2033/PPFGR/outPPFGR_point.dat",
    "WFI2033/NP/outNP_point.dat",
    "WFI2033/NPR/outNPR_point.dat",
    "WFI2033/NPG/outNPG_point.dat",
    "WFI2033/NPGR/outNPGR_point.dat",
    "WFI2033/NPF/outNPF_point.dat",
    "WFI2033/NPFR/outNPFR_point.dat",
    "WFI2033/NPFG/outNPFG_point.dat",
    "WFI2033/NPFGR/outNPFGR_point.dat"
]

Fits_files = [
    "WFI2033/SP/outSP_lens.fits",
    "WFI2033/SPR/outSPR_lens.fits",
    "WFI2033/SPG/outSPG_lens.fits",
    "WFI2033/SPGR/outSPGR_lens.fits",
    "WFI2033/SPF/outSPF_lens.fits",
    "WFI2033/SPFR/outSPFR_lens.fits",
    "WFI2033/SPFG/outSPFG_lens.fits",
    "WFI2033/SPFGR/outSPFGR_lens.fits",
    "WFI2033/PP/outPP_lens.fits",
    "WFI2033/PPR/outPPR_lens.fits",
    "WFI2033/PPG/outPPG_lens.fits",
    "WFI2033/PPGR/outPPGR_lens.fits",
    "WFI2033/PPF/outPPF_lens.fits",
    "WFI2033/PPFR/outPPFR_lens.fits",
    "WFI2033/PPFG/outPPFG_lens.fits",
    "WFI2033/PPFGR/outPPFGR_lens.fits",
    "WFI2033/NP/outNP_lens.fits",
    "WFI2033/NPR/outNPR_lens.fits",
    "WFI2033/NPG/outNPG_lens.fits",
    "WFI2033/NPGR/outNPGR_lens.fits",
    "WFI2033/NPF/outNPF_lens.fits",
    "WFI2033/NPFR/outNPFR_lens.fits",
    "WFI2033/NPFG/outNPFG_lens.fits",
    "WFI2033/NPFGR/outNPFGR_lens.fits"
]

Optresult_files = [
    "WFI2033/SP/outSP_optresult.dat",
    "WFI2033/SPR/outSPR_optresult.dat",
    "WFI2033/SPG/outSPG_optresult.dat",
    "WFI2033/SPGR/outSPGR_optresult.dat",
    "WFI2033/SPF/outSPF_optresult.dat",
    "WFI2033/SPFR/outSPFR_optresult.dat",
    "WFI2033/SPFG/outSPFG_optresult.dat",
    "WFI2033/SPFGR/outSPFGR_optresult.dat",
    "WFI2033/PP/outPP_optresult.dat",
    "WFI2033/PPR/outPPR_optresult.dat",
    "WFI2033/PPG/outPPG_optresult.dat",
    "WFI2033/PPGR/outPPGR_optresult.dat",
    "WFI2033/PPF/outPPF_optresult.dat",
    "WFI2033/PPFR/outPPFR_optresult.dat",
    "WFI2033/PPFG/outPPFG_optresult.dat",
    "WFI2033/PPFGR/outPPFGR_optresult.dat",
    "WFI2033/NP/outNP_optresult.dat",
    "WFI2033/NPR/outNPR_optresult.dat",
    "WFI2033/NPG/outNPG_optresult.dat",
    "WFI2033/NPGR/outNPGR_optresult.dat",
    "WFI2033/NPF/outNPF_optresult.dat",
    "WFI2033/NPFR/outNPFR_optresult.dat",
    "WFI2033/NPFG/outNPFG_optresult.dat",
    "WFI2033/NPFGR/outNPFGR_optresult.dat"
]

n = 0.0025  # Example scale factor
results_df, latex_rows = process_files(Input_file, Output_files, Fits_files, Optresult_files, n)

# Print formatted results
for i, output_file in enumerate(Output_files):
    constraint_type = "PF" if "PF" in Optresult_files[i] else "P"
    print(f"Model: {output_file}")
    print(f"Constraint Type: {constraint_type}")
    print(f"Number of data points used: {results_df.iloc[i]['n_data_points']}")
    print(f"Number of free parameters: {results_df.iloc[i]['n_free_params']}")
    print(f"Delta RMS of Image Positions: {results_df.iloc[i]['delta_rms']:.3f}")
    
    # Print individual position deltas (4 positions)
    print(f"Individual Position Deltas (mas):")
    for j, delta in enumerate(results_df.iloc[i]['individual_position_deltas']):
        print(f"  Image {j+1}: {delta:.3f}")
    
    print(f"Delta RMS of flux ratios at the observed position: {results_df.iloc[i]['predicted_flux_ratios_observed_rms']:.3f}")
    print(f"Delta RMS of flux ratios at the predicted position: {results_df.iloc[i]['predicted_flux_ratios_rms']:.3f}")
    
    # Print absolute flux ratio differences (3 ratios, excluding normalization)
    print(f"Absolute Flux Ratio Differences:")
    flux_diff_indices = [0, 2, 3] if len(results_df.iloc[i]['flux_ratio_differences']) == 3 else range(len(results_df.iloc[i]['flux_ratio_differences']))
    for j, diff in enumerate(results_df.iloc[i]['flux_ratio_differences']):
        img_idx = flux_diff_indices[j] + 1
        print(f"  Image {img_idx}: {diff:.3f}")
    
    print(f"Manual Chi²: {results_df.iloc[i]['manual_chi2']:.3f}")
    if results_df.iloc[i]['glafic_chi2'] is not None:
        print(f"Glafic Chi² (final, reference): {results_df.iloc[i]['glafic_chi2']:.3f}")
        print(f"Chi² Difference (Manual - Glafic): {results_df.iloc[i]['manual_chi2'] - results_df.iloc[i]['glafic_chi2']:.3f}")
    else:
        print(f"Glafic Chi² (reference): Not found")
    print(f"BIC (using manual Chi²): {results_df.iloc[i]['bic']:.3f}")
    
    # Manual verification of BIC calculation
    manual_bic_check = results_df.iloc[i]['manual_chi2'] + results_df.iloc[i]['n_free_params'] * np.log(results_df.iloc[i]['n_data_points'])
    print(f"Manual BIC verification: {results_df.iloc[i]['manual_chi2']:.3f} + {results_df.iloc[i]['n_free_params']} × ln({results_df.iloc[i]['n_data_points']}) = {manual_bic_check:.3f}")
    print()

# Print LaTeX-formatted rows
print("LaTeX Table Rows:")
print("Format: Model & Constraint & Pos1 & Pos2 & Pos3 & Pos4 & Overall ΔRMS & FluxDiff1 & FluxDiff2 & FluxDiff3 & Flux RMS & BIC")
for row in latex_rows:
    print(row)


# In[13]:


# NEW LATEX AUTOMATED CODE WITH MANUAL CHI^2 CALCULATION with individual position delta RMS and flux ratio differences
# SDSSJ1330
import numpy as np
import pandas as pd
from astropy.io import fits
import os

# Function to calculate distance
def calculate_distance(x1, y1, x2, y2):
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

# Function to read observed positions, flux ratios, and uncertainties
def read_observed_positions(input_file):
    observed_positions = []
    observed_flux_ratios = []
    position_uncertainties = []
    flux_uncertainties = []
    with open(input_file, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line
        for line in lines:
            columns = line.split()
            if len(columns) >= 5:
                x1 = float(columns[0])
                y1 = float(columns[1])
                flux_ratio = float(columns[2])
                pos_err = float(columns[3])
                flux_err = float(columns[4])
                observed_positions.append((x1, y1))
                observed_flux_ratios.append(flux_ratio)
                position_uncertainties.append(pos_err)
                flux_uncertainties.append(flux_err)
    return observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties

# Function to read predicted positions and magnifications
def read_predicted_positions(output_file):
    predicted_positions = []
    predicted_magnifications = []
    with open(output_file, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line
        for line in lines:
            if line.startswith("#"):
                continue
            columns = line.split()
            if len(columns) >= 3:
                x2 = float(columns[0])
                y2 = float(columns[1])
                M2 = abs(float(columns[2]))
                predicted_positions.append((x2, y2))
                predicted_magnifications.append(M2)
    return predicted_positions, predicted_magnifications

# Function to pair predicted and observed positions
def pair_positions(observed_positions, predicted_positions, predicted_magnifications):
    paired_positions = []
    for observed_position in observed_positions:
        min_distance = float('inf')
        closest_index = None
        for i, predicted_position in enumerate(predicted_positions):
            distance = calculate_distance(*observed_position, *predicted_position)
            if distance < min_distance:
                min_distance = distance
                closest_index = i
        paired_positions.append((observed_position, predicted_positions[closest_index], predicted_magnifications[closest_index]))
    return paired_positions

# Function to calculate manual chi² with proper uncertainties
def calculate_manual_chi2(observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties, 
                         predicted_positions, predicted_magnifications, green_index=1):
    # Pair positions
    paired_positions = pair_positions(observed_positions, predicted_positions, predicted_magnifications)
    
    # Order predicted positions and magnifications based on observed order
    predicted_positions_ordered = [pair[1] for pair in paired_positions]
    predicted_magnifications_ordered = [pair[2] for pair in paired_positions]
    
    # Calculate predicted flux ratios
    predicted_flux_ratios = [m / predicted_magnifications_ordered[green_index] for m in predicted_magnifications_ordered]
    
    chi2 = 0.0
    
    # Position contributions to chi²
    for i, (observed, predicted, _) in enumerate(paired_positions):
        dx = observed[0] - predicted[0]
        dy = observed[1] - predicted[1]
        pos_err = position_uncertainties[i]
        chi2 += (dx/pos_err)**2 + (dy/pos_err)**2
    
    # Flux ratio contributions to chi² (excluding the normalization image)
    for i, (pred_flux, obs_flux) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)):
        if i != green_index:  # Exclude the green image used for normalization
            flux_err = flux_uncertainties[i]
            chi2 += ((pred_flux - obs_flux)/flux_err)**2
    
    return chi2

# Function to extract the final chi^2 from optresult.dat file
def extract_glafic_chi2(optresult_path):
    chi2_value = None
    with open(optresult_path, 'r') as file:
        lines = file.readlines()
    
    # Find the last occurrence of chi^2 (final/best result)
    for line in lines:
        if line.startswith("chi^2"):
            chi2_value = float(line.split()[2])
    
    return chi2_value

# Function to process file group
def process_file_group(output_file, fits_file, input_file, n):
    observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties = read_observed_positions(input_file)
    predicted_positions, predicted_magnifications = read_predicted_positions(output_file)
    
    paired_positions = pair_positions(observed_positions, predicted_positions, predicted_magnifications)
    
    # Calculate Delta RMS of Image Positions correctly
    squared_distances = 0.0
    individual_distances = []  # Store individual distances for each image
    for observed, predicted, _ in paired_positions:
        distance = calculate_distance(*observed, *predicted)
        individual_distances.append(distance * 1000)  # Convert to mas
        squared_distances += distance ** 2  # Summing squared distances
    
    mean_squared_distances = squared_distances / len(paired_positions)
    delta_rms = np.sqrt(mean_squared_distances)  # Taking square root of the mean of squared distances
    delta_rms *= 1000  # Convert to mas
    
    # Order predicted positions and magnifications based on observed order
    predicted_positions_ordered = [pair[1] for pair in paired_positions]
    predicted_magnifications_ordered = [pair[2] for pair in paired_positions]
    
    # Determine the index of the "Green Image"
    green_index = 0  # CHOOSE WHICH IMAGE TO NORMALISE BY HERE !!!!!!!!!!!!!!!!!!!!!!!

    # Predicted Flux Ratios
    predicted_flux_ratios = [m / predicted_magnifications_ordered[green_index] for m in predicted_magnifications_ordered]

    # Calculate absolute differences between predicted and observed flux ratios (excluding normalization image)
    flux_ratio_differences = []
    for i, (pred_flux, obs_flux) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)):
        if i != green_index:  # Exclude the green image used for normalization
            flux_ratio_differences.append(abs(pred_flux - obs_flux))

    # Calculate manual chi²
    manual_chi2 = calculate_manual_chi2(observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties,
                                       predicted_positions, predicted_magnifications, green_index)

    with fits.open(fits_file) as hdul:
        mag = hdul[0].data[6]
    
    pixel_coords = [(int((x + mag.shape[1] // 2 * n) / n), int((y + mag.shape[0] // 2 * n) / n)) for x, y in observed_positions]
    values = [mag[y, x] for x, y in pixel_coords]
    
    abs_values = [1 / abs(value) for value in values]
    predicted_flux_ratios_observed = [abs_values[i] / abs_values[green_index] for i in range(len(abs_values))]

    # Calculate RMS excluding the green image
    squared_diffs_pred = [(pred - obs) ** 2 for i, (pred, obs) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)) if i != green_index]
    predicted_flux_ratios_rms = np.sqrt(sum(squared_diffs_pred) / len(squared_diffs_pred))

    squared_diffs_obs = [(pred - obs) ** 2 for i, (pred, obs) in enumerate(zip(predicted_flux_ratios_observed, observed_flux_ratios)) if i != green_index]
    predicted_flux_ratios_observed_rms = np.sqrt(sum(squared_diffs_obs) / len(squared_diffs_obs))

    return {
        'delta_rms': delta_rms,
        'individual_position_deltas': individual_distances,  # Delta RMS for each of the 4 positions
        'flux_ratio_differences': flux_ratio_differences,    # Absolute differences for 3 flux ratios
        'predicted_flux_ratios_rms': predicted_flux_ratios_rms,
        'predicted_flux_ratios_observed_rms': predicted_flux_ratios_observed_rms,
        'manual_chi2': manual_chi2
    }

# Function to extract chi^2 and lens models from optresult.dat file (following Code 1's approach exactly)
def extract_chi2_and_lens_profiles(optresult_path):
    chi2_value = None
    profiles_used = []
    sie_count = 0

    with open(optresult_path, 'r') as file:
        lines = file.readlines()

    found_chi2 = False
    for line in lines:
        if line.startswith("chi^2"):
            chi2_value = float(line.split()[2])
            found_chi2 = True
            profiles_used = []
            sie_count = 0
        elif found_chi2:
            if line.startswith("lens   sie"):
                profiles_used.append("SIE")
                sie_count += 1
            elif line.startswith("lens   pow"):
                profiles_used.append("POW")
            elif line.startswith("lens   anfw"):
                profiles_used.append("NFW")
            elif line.startswith("lens   pert"):
                profiles_used.append("Shear")

    if chi2_value is None:
        raise ValueError("No chi^2 value found in the file.")

    profiles_used = list(set(profiles_used))
    return chi2_value, profiles_used, sie_count

# Function to calculate free parameters based on lens profiles (identical to Code 1)
def calculate_free_parameters(profiles_used, sie_count):
    profile_params = {
        "SIE": 5,
        "POW": 6,
        "NFW": 6,
        "Shear": 5
    }
    total_params = sum(profile_params[profile] for profile in profiles_used)
    if sie_count == 2:
        total_params += 5
    return total_params

# Function to calculate BIC (Bayesian Information Criterion)
def calculate_bic(chi2, n_data_points, n_free_params):
    return chi2 + n_free_params * np.log(n_data_points)

# Function to compute BIC using manual chi² but following Code 1's approach for free parameters
def compute_bic_with_tracking(manual_chi2, optresult_path):
    # Use Code 1's approach: extract profiles from optresult.dat file
    glafic_chi2, profiles_used, sie_count = extract_chi2_and_lens_profiles(optresult_path)
    n_free_params = calculate_free_parameters(profiles_used, sie_count)
    
    # Properly adjust n_data_points based on constraint type
    if "PF" in optresult_path:
        n_data_points = 4 * 2 + 3  # positions + flux ratios
    else:
        n_data_points = 4 * 2  # positions only
    
    # Always use manual chi² for BIC calculation
    bic_value = calculate_bic(manual_chi2, n_data_points, n_free_params)
    
    return bic_value, glafic_chi2, n_data_points, n_free_params

# Main function to process multiple files
def process_files(input_file, output_files, fits_files, optresult_files, n):
    results = []
    latex_rows = []
    bics = []  # Store BIC values separately
    glafic_chi2_values = []  # Store glafic chi² for reference
    data_points_used = []  # Store number of data points used
    free_params_used = []  # Store number of free parameters used

    for output_file, fits_file, optresult_file in zip(output_files, fits_files, optresult_files):
        # Process output and fits files
        result = process_file_group(output_file, fits_file, input_file, n)
        results.append(result)
        
        # Calculate BIC using manual chi² but track both values
        bic, glafic_chi2, n_data_points, n_free_params = compute_bic_with_tracking(result['manual_chi2'], optresult_file)
        bics.append(bic)
        glafic_chi2_values.append(glafic_chi2)
        data_points_used.append(n_data_points)
        free_params_used.append(n_free_params)
        
        # Create LaTeX row with additional information
        model_name, constraint_name = determine_model_and_constraint(output_file)
        # Format individual position deltas and flux ratio differences for LaTeX
        pos_deltas_str = " & ".join([f"{delta:.3f}" for delta in result['individual_position_deltas']])
        flux_diffs_str = " & ".join([f"{diff:.3f}" for diff in result['flux_ratio_differences']])
        
        latex_row = f"{model_name} & {constraint_name} & {pos_deltas_str} & {result['delta_rms']:.1f} & {flux_diffs_str} & {result['predicted_flux_ratios_rms']:.2f} & {bic:.1f} \\\\"
        latex_rows.append(latex_row)
    
    # Combine results into a DataFrame
    df_results = pd.DataFrame(results)
    df_results['bic'] = bics  # Add BIC values to the DataFrame
    df_results['glafic_chi2'] = glafic_chi2_values  # Add glafic chi² for reference
    df_results['n_data_points'] = data_points_used  # Add number of data points used
    df_results['n_free_params'] = free_params_used  # Add number of free parameters used
    
    return df_results, latex_rows

# Determine model and constraint names from the filename
def determine_model_and_constraint(file_name):
    model_map = {'S': 'SIE', 'P': 'POW', 'N': 'NFW'}
    
    base_name = os.path.basename(file_name)
    model_key = base_name[3]
    constraint_key = base_name[4:6] if base_name[4:6] in ['PF'] else base_name[4]
    
    model_name = model_map.get(model_key, '')
    if 'R' in base_name:
        model_name += ' + Shear'
    if 'G' in base_name:
        model_name += ' + G2'
    
    constraint_name = 'QSO Pos + FR' if constraint_key == 'PF' else 'QSO Pos'
    
    return model_name, constraint_name

# Example usage with filenames
Input_file = "SDSSJ1330/SPFC/Robs_point.dat"

Output_files = [
    "SDSSJ1330/SPC/outSP_point.dat",
    "SDSSJ1330/SPC/outSPR_point.dat",
    "SDSSJ1330/SPFC/outSPF_point.dat",
    "SDSSJ1330/SPFC/outSPFR_point.dat",
    "SDSSJ1330/PPC/outPP_point.dat",
    "SDSSJ1330/PPC/outPPR_point.dat",
    "SDSSJ1330/PPFC/outPPF_point.dat",
    "SDSSJ1330/PPFC/outPPFR_point.dat",
    "SDSSJ1330/NPC/outNP_point.dat",
    "SDSSJ1330/NPC/outNPR_point.dat",
    "SDSSJ1330/NPFC/outNPF_point.dat",
    "SDSSJ1330/NPFC/outNPFR_point.dat",
]

Fits_files = [
    "SDSSJ1330/SPC/outSP_lens.fits",
    "SDSSJ1330/SPC/outSPR_lens.fits",
    "SDSSJ1330/SPFC/outSPF_lens.fits",
    "SDSSJ1330/SPFC/outSPFR_lens.fits",
    "SDSSJ1330/PPC/outPP_lens.fits",
    "SDSSJ1330/PPC/outPPR_lens.fits",
    "SDSSJ1330/PPFC/outPPF_lens.fits",
    "SDSSJ1330/PPFC/outPPFR_lens.fits",
    "SDSSJ1330/NPC/outNP_lens.fits",
    "SDSSJ1330/NPC/outNPR_lens.fits",
    "SDSSJ1330/NPFC/outNPF_lens.fits",
    "SDSSJ1330/NPFC/outNPFR_lens.fits",
]

Optresult_files = [
    "SDSSJ1330/SPC/outSP_optresult.dat",
    "SDSSJ1330/SPC/outSPR_optresult.dat",
    "SDSSJ1330/SPFC/outSPF_optresult.dat",
    "SDSSJ1330/SPFC/outSPFR_optresult.dat",
    "SDSSJ1330/PPC/outPP_optresult.dat",
    "SDSSJ1330/PPC/outPPR_optresult.dat",
    "SDSSJ1330/PPFC/outPPF_optresult.dat",
    "SDSSJ1330/PPFC/outPPFR_optresult.dat",
    "SDSSJ1330/NPC/outNP_optresult.dat",
    "SDSSJ1330/NPC/outNPR_optresult.dat",
    "SDSSJ1330/NPFC/outNPF_optresult.dat",
    "SDSSJ1330/NPFC/outNPFR_optresult.dat",
]

n = 0.0025  # Example scale factor
results_df, latex_rows = process_files(Input_file, Output_files, Fits_files, Optresult_files, n)

# Print formatted results
for i, output_file in enumerate(Output_files):
    constraint_type = "PF" if "PF" in Optresult_files[i] else "P"
    print(f"Model: {output_file}")
    print(f"Constraint Type: {constraint_type}")
    print(f"Number of data points used: {results_df.iloc[i]['n_data_points']}")
    print(f"Number of free parameters: {results_df.iloc[i]['n_free_params']}")
    print(f"Delta RMS of Image Positions: {results_df.iloc[i]['delta_rms']:.3f}")
    
    # Print individual position deltas (4 positions)
    print(f"Individual Position Deltas (mas):")
    for j, delta in enumerate(results_df.iloc[i]['individual_position_deltas']):
        print(f"  Image {j+1}: {delta:.3f}")
    
    print(f"Delta RMS of flux ratios at the observed position: {results_df.iloc[i]['predicted_flux_ratios_observed_rms']:.3f}")
    print(f"Delta RMS of flux ratios at the predicted position: {results_df.iloc[i]['predicted_flux_ratios_rms']:.3f}")
    
    # Print absolute flux ratio differences (3 ratios, excluding normalization)
    print(f"Absolute Flux Ratio Differences:")
    flux_diff_indices = [0, 2, 3] if len(results_df.iloc[i]['flux_ratio_differences']) == 3 else range(len(results_df.iloc[i]['flux_ratio_differences']))
    for j, diff in enumerate(results_df.iloc[i]['flux_ratio_differences']):
        img_idx = flux_diff_indices[j] + 1
        print(f"  Image {img_idx}: {diff:.3f}")
    
    print(f"Manual Chi²: {results_df.iloc[i]['manual_chi2']:.3f}")
    if results_df.iloc[i]['glafic_chi2'] is not None:
        print(f"Glafic Chi² (final, reference): {results_df.iloc[i]['glafic_chi2']:.3f}")
        print(f"Chi² Difference (Manual - Glafic): {results_df.iloc[i]['manual_chi2'] - results_df.iloc[i]['glafic_chi2']:.3f}")
    else:
        print(f"Glafic Chi² (reference): Not found")
    print(f"BIC (using manual Chi²): {results_df.iloc[i]['bic']:.3f}")
    
    # Manual verification of BIC calculation
    manual_bic_check = results_df.iloc[i]['manual_chi2'] + results_df.iloc[i]['n_free_params'] * np.log(results_df.iloc[i]['n_data_points'])
    print(f"Manual BIC verification: {results_df.iloc[i]['manual_chi2']:.3f} + {results_df.iloc[i]['n_free_params']} × ln({results_df.iloc[i]['n_data_points']}) = {manual_bic_check:.3f}")
    print()

# Print LaTeX-formatted rows
print("LaTeX Table Rows:")
print("Format: Model & Constraint & Pos1 & Pos2 & Pos3 & Pos4 & Overall ΔRMS & FluxDiff1 & FluxDiff2 & FluxDiff3 & Flux RMS & BIC")
for row in latex_rows:
    print(row)


# In[7]:


# NEW LATEX AUTOMATED CODE WITH MANUAL CHI^2 CALCULATION with individual position delta RMS and flux ratio differences
# WFI2026
import numpy as np
import pandas as pd
from astropy.io import fits
import os

# Function to calculate distance
def calculate_distance(x1, y1, x2, y2):
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

# Function to read observed positions, flux ratios, and uncertainties
def read_observed_positions(input_file):
    observed_positions = []
    observed_flux_ratios = []
    position_uncertainties = []
    flux_uncertainties = []
    with open(input_file, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line
        for line in lines:
            columns = line.split()
            if len(columns) >= 5:
                x1 = float(columns[0])
                y1 = float(columns[1])
                flux_ratio = float(columns[2])
                pos_err = float(columns[3])
                flux_err = float(columns[4])
                observed_positions.append((x1, y1))
                observed_flux_ratios.append(flux_ratio)
                position_uncertainties.append(pos_err)
                flux_uncertainties.append(flux_err)
    return observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties

# Function to read predicted positions and magnifications
def read_predicted_positions(output_file):
    predicted_positions = []
    predicted_magnifications = []
    with open(output_file, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line
        for line in lines:
            if line.startswith("#"):
                continue
            columns = line.split()
            if len(columns) >= 3:
                x2 = float(columns[0])
                y2 = float(columns[1])
                M2 = abs(float(columns[2]))
                predicted_positions.append((x2, y2))
                predicted_magnifications.append(M2)
    return predicted_positions, predicted_magnifications

# Function to pair predicted and observed positions
def pair_positions(observed_positions, predicted_positions, predicted_magnifications):
    paired_positions = []
    for observed_position in observed_positions:
        min_distance = float('inf')
        closest_index = None
        for i, predicted_position in enumerate(predicted_positions):
            distance = calculate_distance(*observed_position, *predicted_position)
            if distance < min_distance:
                min_distance = distance
                closest_index = i
        paired_positions.append((observed_position, predicted_positions[closest_index], predicted_magnifications[closest_index]))
    return paired_positions

# Function to calculate manual chi² with proper uncertainties
def calculate_manual_chi2(observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties, 
                         predicted_positions, predicted_magnifications, green_index=1):
    # Pair positions
    paired_positions = pair_positions(observed_positions, predicted_positions, predicted_magnifications)
    
    # Order predicted positions and magnifications based on observed order
    predicted_positions_ordered = [pair[1] for pair in paired_positions]
    predicted_magnifications_ordered = [pair[2] for pair in paired_positions]
    
    # Calculate predicted flux ratios
    predicted_flux_ratios = [m / predicted_magnifications_ordered[green_index] for m in predicted_magnifications_ordered]
    
    chi2 = 0.0
    
    # Position contributions to chi²
    for i, (observed, predicted, _) in enumerate(paired_positions):
        dx = observed[0] - predicted[0]
        dy = observed[1] - predicted[1]
        pos_err = position_uncertainties[i]
        chi2 += (dx/pos_err)**2 + (dy/pos_err)**2
    
    # Flux ratio contributions to chi² (excluding the normalization image)
    for i, (pred_flux, obs_flux) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)):
        if i != green_index:  # Exclude the green image used for normalization
            flux_err = flux_uncertainties[i]
            chi2 += ((pred_flux - obs_flux)/flux_err)**2
    
    return chi2

# Function to extract the final chi^2 from optresult.dat file
def extract_glafic_chi2(optresult_path):
    chi2_value = None
    with open(optresult_path, 'r') as file:
        lines = file.readlines()
    
    # Find the last occurrence of chi^2 (final/best result)
    for line in lines:
        if line.startswith("chi^2"):
            chi2_value = float(line.split()[2])
    
    return chi2_value

# Function to process file group
def process_file_group(output_file, fits_file, input_file, n):
    observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties = read_observed_positions(input_file)
    predicted_positions, predicted_magnifications = read_predicted_positions(output_file)
    
    paired_positions = pair_positions(observed_positions, predicted_positions, predicted_magnifications)
    
    # Calculate Delta RMS of Image Positions correctly
    squared_distances = 0.0
    individual_distances = []  # Store individual distances for each image
    for observed, predicted, _ in paired_positions:
        distance = calculate_distance(*observed, *predicted)
        individual_distances.append(distance * 1000)  # Convert to mas
        squared_distances += distance ** 2  # Summing squared distances
    
    mean_squared_distances = squared_distances / len(paired_positions)
    delta_rms = np.sqrt(mean_squared_distances)  # Taking square root of the mean of squared distances
    delta_rms *= 1000  # Convert to mas
    
    # Order predicted positions and magnifications based on observed order
    predicted_positions_ordered = [pair[1] for pair in paired_positions]
    predicted_magnifications_ordered = [pair[2] for pair in paired_positions]
    
    # Determine the index of the "Green Image"
    green_index = 0  # CHOOSE WHICH IMAGE TO NORMALISE BY HERE !!!!!!!!!!!!!!!!!!!!!!!

    # Predicted Flux Ratios
    predicted_flux_ratios = [m / predicted_magnifications_ordered[green_index] for m in predicted_magnifications_ordered]

    # Calculate absolute differences between predicted and observed flux ratios (excluding normalization image)
    flux_ratio_differences = []
    for i, (pred_flux, obs_flux) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)):
        if i != green_index:  # Exclude the green image used for normalization
            flux_ratio_differences.append(abs(pred_flux - obs_flux))

    # Calculate manual chi²
    manual_chi2 = calculate_manual_chi2(observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties,
                                       predicted_positions, predicted_magnifications, green_index)

    with fits.open(fits_file) as hdul:
        mag = hdul[0].data[6]
    
    pixel_coords = [(int((x + mag.shape[1] // 2 * n) / n), int((y + mag.shape[0] // 2 * n) / n)) for x, y in observed_positions]
    values = [mag[y, x] for x, y in pixel_coords]
    
    abs_values = [1 / abs(value) for value in values]
    predicted_flux_ratios_observed = [abs_values[i] / abs_values[green_index] for i in range(len(abs_values))]

    # Calculate RMS excluding the green image
    squared_diffs_pred = [(pred - obs) ** 2 for i, (pred, obs) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)) if i != green_index]
    predicted_flux_ratios_rms = np.sqrt(sum(squared_diffs_pred) / len(squared_diffs_pred))

    squared_diffs_obs = [(pred - obs) ** 2 for i, (pred, obs) in enumerate(zip(predicted_flux_ratios_observed, observed_flux_ratios)) if i != green_index]
    predicted_flux_ratios_observed_rms = np.sqrt(sum(squared_diffs_obs) / len(squared_diffs_obs))

    return {
        'delta_rms': delta_rms,
        'individual_position_deltas': individual_distances,  # Delta RMS for each of the 4 positions
        'flux_ratio_differences': flux_ratio_differences,    # Absolute differences for 3 flux ratios
        'predicted_flux_ratios_rms': predicted_flux_ratios_rms,
        'predicted_flux_ratios_observed_rms': predicted_flux_ratios_observed_rms,
        'manual_chi2': manual_chi2
    }

# Function to extract chi^2 and lens models from optresult.dat file (following Code 1's approach exactly)
def extract_chi2_and_lens_profiles(optresult_path):
    chi2_value = None
    profiles_used = []
    sie_count = 0

    with open(optresult_path, 'r') as file:
        lines = file.readlines()

    found_chi2 = False
    for line in lines:
        if line.startswith("chi^2"):
            chi2_value = float(line.split()[2])
            found_chi2 = True
            profiles_used = []
            sie_count = 0
        elif found_chi2:
            if line.startswith("lens   sie"):
                profiles_used.append("SIE")
                sie_count += 1
            elif line.startswith("lens   pow"):
                profiles_used.append("POW")
            elif line.startswith("lens   anfw"):
                profiles_used.append("NFW")
            elif line.startswith("lens   pert"):
                profiles_used.append("Shear")

    if chi2_value is None:
        raise ValueError("No chi^2 value found in the file.")

    profiles_used = list(set(profiles_used))
    return chi2_value, profiles_used, sie_count

# Function to calculate free parameters based on lens profiles (identical to Code 1)
def calculate_free_parameters(profiles_used, sie_count):
    profile_params = {
        "SIE": 5,
        "POW": 6,
        "NFW": 6,
        "Shear": 5
    }
    total_params = sum(profile_params[profile] for profile in profiles_used)
    if sie_count == 2:
        total_params += 5
    return total_params

# Function to calculate BIC (Bayesian Information Criterion)
def calculate_bic(chi2, n_data_points, n_free_params):
    return chi2 + n_free_params * np.log(n_data_points)

# Function to compute BIC using manual chi² but following Code 1's approach for free parameters
def compute_bic_with_tracking(manual_chi2, optresult_path):
    # Use Code 1's approach: extract profiles from optresult.dat file
    glafic_chi2, profiles_used, sie_count = extract_chi2_and_lens_profiles(optresult_path)
    n_free_params = calculate_free_parameters(profiles_used, sie_count)
    
    # Properly adjust n_data_points based on constraint type
    if "PF" in optresult_path:
        n_data_points = 4 * 2 + 3  # positions + flux ratios
    else:
        n_data_points = 4 * 2  # positions only
    
    # Always use manual chi² for BIC calculation
    bic_value = calculate_bic(manual_chi2, n_data_points, n_free_params)
    
    return bic_value, glafic_chi2, n_data_points, n_free_params

# Main function to process multiple files
def process_files(input_file, output_files, fits_files, optresult_files, n):
    results = []
    latex_rows = []
    bics = []  # Store BIC values separately
    glafic_chi2_values = []  # Store glafic chi² for reference
    data_points_used = []  # Store number of data points used
    free_params_used = []  # Store number of free parameters used

    for output_file, fits_file, optresult_file in zip(output_files, fits_files, optresult_files):
        # Process output and fits files
        result = process_file_group(output_file, fits_file, input_file, n)
        results.append(result)
        
        # Calculate BIC using manual chi² but track both values
        bic, glafic_chi2, n_data_points, n_free_params = compute_bic_with_tracking(result['manual_chi2'], optresult_file)
        bics.append(bic)
        glafic_chi2_values.append(glafic_chi2)
        data_points_used.append(n_data_points)
        free_params_used.append(n_free_params)
        
        # Create LaTeX row with additional information
        model_name, constraint_name = determine_model_and_constraint(output_file)
        # Format individual position deltas and flux ratio differences for LaTeX
        pos_deltas_str = " & ".join([f"{delta:.3f}" for delta in result['individual_position_deltas']])
        flux_diffs_str = " & ".join([f"{diff:.3f}" for diff in result['flux_ratio_differences']])
        
        latex_row = f"{model_name} & {constraint_name} & {pos_deltas_str} & {result['delta_rms']:.1f} & {flux_diffs_str} & {result['predicted_flux_ratios_rms']:.2f} & {bic:.1f} \\\\"
        latex_rows.append(latex_row)
    
    # Combine results into a DataFrame
    df_results = pd.DataFrame(results)
    df_results['bic'] = bics  # Add BIC values to the DataFrame
    df_results['glafic_chi2'] = glafic_chi2_values  # Add glafic chi² for reference
    df_results['n_data_points'] = data_points_used  # Add number of data points used
    df_results['n_free_params'] = free_params_used  # Add number of free parameters used
    
    return df_results, latex_rows

# Determine model and constraint names from the filename
def determine_model_and_constraint(file_name):
    model_map = {'S': 'SIE', 'P': 'POW', 'N': 'NFW'}
    
    base_name = os.path.basename(file_name)
    model_key = base_name[3]
    constraint_key = base_name[4:6] if base_name[4:6] in ['PF'] else base_name[4]
    
    model_name = model_map.get(model_key, '')
    if 'R' in base_name:
        model_name += ' + Shear'
    if 'G' in base_name:
        model_name += ' + G2'
    
    constraint_name = 'QSO Pos + FR' if constraint_key == 'PF' else 'QSO Pos'
    
    return model_name, constraint_name

# Example usage with filenames
Input_file = "WFI2026/SPFC/Wobs_point.dat"

Output_files = [
    "WFI2026/SPC/outSP_point.dat",
    "WFI2026/SPC/outSPR_point.dat",
    "WFI2026/SPFC/outSPF_point.dat",
    "WFI2026/SPFC/outSPFR_point.dat",
    "WFI2026/PPC/outPP_point.dat",
    "WFI2026/PPC/outPPR_point.dat",
    "WFI2026/PPFC/outPPF_point.dat",
    "WFI2026/PPFC/outPPFR_point.dat",
    "WFI2026/NPC/outNP_point.dat",
    "WFI2026/NPC/outNPR_point.dat",
    "WFI2026/NPFC/outNPF_point.dat",
    "WFI2026/NPFC/outNPFR_point.dat",
]

Fits_files = [
    "WFI2026/SPC/outSP_lens.fits",
    "WFI2026/SPC/outSPR_lens.fits",
    "WFI2026/SPFC/outSPF_lens.fits",
    "WFI2026/SPFC/outSPFR_lens.fits",
    "WFI2026/PPC/outPP_lens.fits",
    "WFI2026/PPC/outPPR_lens.fits",
    "WFI2026/PPFC/outPPF_lens.fits",
    "WFI2026/PPFC/outPPFR_lens.fits",
    "WFI2026/NPC/outNP_lens.fits",
    "WFI2026/NPC/outNPR_lens.fits",
    "WFI2026/NPFC/outNPF_lens.fits",
    "WFI2026/NPFC/outNPFR_lens.fits",
]

Optresult_files = [
    "WFI2026/SPC/outSP_optresult.dat",
    "WFI2026/SPC/outSPR_optresult.dat",
    "WFI2026/SPFC/outSPF_optresult.dat",
    "WFI2026/SPFC/outSPFR_optresult.dat",
    "WFI2026/PPC/outPP_optresult.dat",
    "WFI2026/PPC/outPPR_optresult.dat",
    "WFI2026/PPFC/outPPF_optresult.dat",
    "WFI2026/PPFC/outPPFR_optresult.dat",
    "WFI2026/NPC/outNP_optresult.dat",
    "WFI2026/NPC/outNPR_optresult.dat",
    "WFI2026/NPFC/outNPF_optresult.dat",
    "WFI2026/NPFC/outNPFR_optresult.dat",
]

n = 0.0025  # Example scale factor
results_df, latex_rows = process_files(Input_file, Output_files, Fits_files, Optresult_files, n)

# Print formatted results
for i, output_file in enumerate(Output_files):
    constraint_type = "PF" if "PF" in Optresult_files[i] else "P"
    print(f"Model: {output_file}")
    print(f"Constraint Type: {constraint_type}")
    print(f"Number of data points used: {results_df.iloc[i]['n_data_points']}")
    print(f"Number of free parameters: {results_df.iloc[i]['n_free_params']}")
    print(f"Delta RMS of Image Positions: {results_df.iloc[i]['delta_rms']:.3f}")
    
    # Print individual position deltas (4 positions)
    print(f"Individual Position Deltas (mas):")
    for j, delta in enumerate(results_df.iloc[i]['individual_position_deltas']):
        print(f"  Image {j+1}: {delta:.3f}")
    
    print(f"Delta RMS of flux ratios at the observed position: {results_df.iloc[i]['predicted_flux_ratios_observed_rms']:.3f}")
    print(f"Delta RMS of flux ratios at the predicted position: {results_df.iloc[i]['predicted_flux_ratios_rms']:.3f}")
    
    # Print absolute flux ratio differences (3 ratios, excluding normalization)
    print(f"Absolute Flux Ratio Differences:")
    flux_diff_indices = [0, 2, 3] if len(results_df.iloc[i]['flux_ratio_differences']) == 3 else range(len(results_df.iloc[i]['flux_ratio_differences']))
    for j, diff in enumerate(results_df.iloc[i]['flux_ratio_differences']):
        img_idx = flux_diff_indices[j] + 1
        print(f"  Image {img_idx}: {diff:.3f}")
    
    print(f"Manual Chi²: {results_df.iloc[i]['manual_chi2']:.3f}")
    if results_df.iloc[i]['glafic_chi2'] is not None:
        print(f"Glafic Chi² (final, reference): {results_df.iloc[i]['glafic_chi2']:.3f}")
        print(f"Chi² Difference (Manual - Glafic): {results_df.iloc[i]['manual_chi2'] - results_df.iloc[i]['glafic_chi2']:.3f}")
    else:
        print(f"Glafic Chi² (reference): Not found")
    print(f"BIC (using manual Chi²): {results_df.iloc[i]['bic']:.3f}")
    
    # Manual verification of BIC calculation
    manual_bic_check = results_df.iloc[i]['manual_chi2'] + results_df.iloc[i]['n_free_params'] * np.log(results_df.iloc[i]['n_data_points'])
    print(f"Manual BIC verification: {results_df.iloc[i]['manual_chi2']:.3f} + {results_df.iloc[i]['n_free_params']} × ln({results_df.iloc[i]['n_data_points']}) = {manual_bic_check:.3f}")
    print()

# Print LaTeX-formatted rows
print("LaTeX Table Rows:")
print("Format: Model & Constraint & Pos1 & Pos2 & Pos3 & Pos4 & Overall ΔRMS & FluxDiff1 & FluxDiff2 & FluxDiff3 & Flux RMS & BIC")
for row in latex_rows:
    print(row)


# In[1]:


# NEW LATEX AUTOMATED CODE WITH MANUAL CHI^2 CALCULATION with individual position delta RMS and flux ratio differences
# WFI2026
import numpy as np
import pandas as pd
from astropy.io import fits
import os

# Function to calculate distance
def calculate_distance(x1, y1, x2, y2):
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

# Function to read observed positions, flux ratios, and uncertainties
def read_observed_positions(input_file):
    observed_positions = []
    observed_flux_ratios = []
    position_uncertainties = []
    flux_uncertainties = []
    with open(input_file, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line
        for line in lines:
            columns = line.split()
            if len(columns) >= 5:
                x1 = float(columns[0])
                y1 = float(columns[1])
                flux_ratio = float(columns[2])
                pos_err = float(columns[3])
                flux_err = float(columns[4])
                observed_positions.append((x1, y1))
                observed_flux_ratios.append(flux_ratio)
                position_uncertainties.append(pos_err)
                flux_uncertainties.append(flux_err)
    return observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties

# Function to read predicted positions and magnifications
def read_predicted_positions(output_file):
    predicted_positions = []
    predicted_magnifications = []
    with open(output_file, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line
        for line in lines:
            if line.startswith("#"):
                continue
            columns = line.split()
            if len(columns) >= 3:
                x2 = float(columns[0])
                y2 = float(columns[1])
                M2 = abs(float(columns[2]))
                predicted_positions.append((x2, y2))
                predicted_magnifications.append(M2)
    return predicted_positions, predicted_magnifications

# Function to pair predicted and observed positions
def pair_positions(observed_positions, predicted_positions, predicted_magnifications):
    paired_positions = []
    for observed_position in observed_positions:
        min_distance = float('inf')
        closest_index = None
        for i, predicted_position in enumerate(predicted_positions):
            distance = calculate_distance(*observed_position, *predicted_position)
            if distance < min_distance:
                min_distance = distance
                closest_index = i
        paired_positions.append((observed_position, predicted_positions[closest_index], predicted_magnifications[closest_index]))
    return paired_positions

# Function to calculate manual chi² with proper uncertainties
def calculate_manual_chi2(observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties, 
                         predicted_positions, predicted_magnifications, green_index=1):
    # Pair positions
    paired_positions = pair_positions(observed_positions, predicted_positions, predicted_magnifications)
    
    # Order predicted positions and magnifications based on observed order
    predicted_positions_ordered = [pair[1] for pair in paired_positions]
    predicted_magnifications_ordered = [pair[2] for pair in paired_positions]
    
    # Calculate predicted flux ratios
    predicted_flux_ratios = [m / predicted_magnifications_ordered[green_index] for m in predicted_magnifications_ordered]
    
    chi2 = 0.0
    
    # Position contributions to chi²
    for i, (observed, predicted, _) in enumerate(paired_positions):
        dx = observed[0] - predicted[0]
        dy = observed[1] - predicted[1]
        pos_err = position_uncertainties[i]
        chi2 += (dx/pos_err)**2 + (dy/pos_err)**2
    
    # Flux ratio contributions to chi² (excluding the normalization image)
    for i, (pred_flux, obs_flux) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)):
        if i != green_index:  # Exclude the green image used for normalization
            flux_err = flux_uncertainties[i]
            chi2 += ((pred_flux - obs_flux)/flux_err)**2
    
    return chi2

# Function to extract the final chi^2 from optresult.dat file
def extract_glafic_chi2(optresult_path):
    chi2_value = None
    with open(optresult_path, 'r') as file:
        lines = file.readlines()
    
    # Find the last occurrence of chi^2 (final/best result)
    for line in lines:
        if line.startswith("chi^2"):
            chi2_value = float(line.split()[2])
    
    return chi2_value

# Function to process file group
def process_file_group(output_file, fits_file, input_file, n):
    observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties = read_observed_positions(input_file)
    predicted_positions, predicted_magnifications = read_predicted_positions(output_file)
    
    paired_positions = pair_positions(observed_positions, predicted_positions, predicted_magnifications)
    
    # Calculate Delta RMS of Image Positions correctly
    squared_distances = 0.0
    individual_distances = []  # Store individual distances for each image
    for observed, predicted, _ in paired_positions:
        distance = calculate_distance(*observed, *predicted)
        individual_distances.append(distance * 1000)  # Convert to mas
        squared_distances += distance ** 2  # Summing squared distances
    
    mean_squared_distances = squared_distances / len(paired_positions)
    delta_rms = np.sqrt(mean_squared_distances)  # Taking square root of the mean of squared distances
    delta_rms *= 1000  # Convert to mas
    
    # Order predicted positions and magnifications based on observed order
    predicted_positions_ordered = [pair[1] for pair in paired_positions]
    predicted_magnifications_ordered = [pair[2] for pair in paired_positions]
    
    # Determine the index of the "Green Image"
    green_index = 0  # CHOOSE WHICH IMAGE TO NORMALISE BY HERE !!!!!!!!!!!!!!!!!!!!!!!

    # Predicted Flux Ratios
    predicted_flux_ratios = [m / predicted_magnifications_ordered[green_index] for m in predicted_magnifications_ordered]

    # Calculate absolute differences between predicted and observed flux ratios (excluding normalization image)
    flux_ratio_differences = []
    for i, (pred_flux, obs_flux) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)):
        if i != green_index:  # Exclude the green image used for normalization
            flux_ratio_differences.append(abs(pred_flux - obs_flux))

    # Calculate manual chi²
    manual_chi2 = calculate_manual_chi2(observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties,
                                       predicted_positions, predicted_magnifications, green_index)

    with fits.open(fits_file) as hdul:
        mag = hdul[0].data[6]
    
    pixel_coords = [(int((x + mag.shape[1] // 2 * n) / n), int((y + mag.shape[0] // 2 * n) / n)) for x, y in observed_positions]
    values = [mag[y, x] for x, y in pixel_coords]
    
    abs_values = [1 / abs(value) for value in values]
    predicted_flux_ratios_observed = [abs_values[i] / abs_values[green_index] for i in range(len(abs_values))]

    # Calculate RMS excluding the green image
    squared_diffs_pred = [(pred - obs) ** 2 for i, (pred, obs) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)) if i != green_index]
    predicted_flux_ratios_rms = np.sqrt(sum(squared_diffs_pred) / len(squared_diffs_pred))

    squared_diffs_obs = [(pred - obs) ** 2 for i, (pred, obs) in enumerate(zip(predicted_flux_ratios_observed, observed_flux_ratios)) if i != green_index]
    predicted_flux_ratios_observed_rms = np.sqrt(sum(squared_diffs_obs) / len(squared_diffs_obs))

    return {
        'delta_rms': delta_rms,
        'individual_position_deltas': individual_distances,  # Delta RMS for each of the 4 positions
        'flux_ratio_differences': flux_ratio_differences,    # Absolute differences for 3 flux ratios
        'predicted_flux_ratios_rms': predicted_flux_ratios_rms,
        'predicted_flux_ratios_observed_rms': predicted_flux_ratios_observed_rms,
        'manual_chi2': manual_chi2
    }

# Function to extract chi^2 and lens models from optresult.dat file (following Code 1's approach exactly)
def extract_chi2_and_lens_profiles(optresult_path):
    chi2_value = None
    profiles_used = []
    sie_count = 0

    with open(optresult_path, 'r') as file:
        lines = file.readlines()

    found_chi2 = False
    for line in lines:
        if line.startswith("chi^2"):
            chi2_value = float(line.split()[2])
            found_chi2 = True
            profiles_used = []
            sie_count = 0
        elif found_chi2:
            if line.startswith("lens   sie"):
                profiles_used.append("SIE")
                sie_count += 1
            elif line.startswith("lens   pow"):
                profiles_used.append("POW")
            elif line.startswith("lens   anfw"):
                profiles_used.append("NFW")
            elif line.startswith("lens   pert"):
                profiles_used.append("Shear")

    if chi2_value is None:
        raise ValueError("No chi^2 value found in the file.")

    profiles_used = list(set(profiles_used))
    return chi2_value, profiles_used, sie_count

# Function to calculate free parameters based on lens profiles (identical to Code 1)
def calculate_free_parameters(profiles_used, sie_count):
    profile_params = {
        "SIE": 5,
        "POW": 6,
        "NFW": 6,
        "Shear": 5
    }
    total_params = sum(profile_params[profile] for profile in profiles_used)
    if sie_count == 2:
        total_params += 5
    return total_params

# Function to calculate BIC (Bayesian Information Criterion)
def calculate_bic(chi2, n_data_points, n_free_params):
    return chi2 + n_free_params * np.log(n_data_points)

# Function to compute BIC using manual chi² but following Code 1's approach for free parameters
def compute_bic_with_tracking(manual_chi2, optresult_path):
    # Use Code 1's approach: extract profiles from optresult.dat file
    glafic_chi2, profiles_used, sie_count = extract_chi2_and_lens_profiles(optresult_path)
    n_free_params = calculate_free_parameters(profiles_used, sie_count)
    
    # Properly adjust n_data_points based on constraint type
    if "PF" in optresult_path:
        n_data_points = 4 * 2 + 3  # positions + flux ratios
    else:
        n_data_points = 4 * 2  # positions only
    
    # Always use manual chi² for BIC calculation
    bic_value = calculate_bic(manual_chi2, n_data_points, n_free_params)
    
    return bic_value, glafic_chi2, n_data_points, n_free_params

# Main function to process multiple files
def process_files(input_file, output_files, fits_files, optresult_files, n):
    results = []
    latex_rows = []
    bics = []  # Store BIC values separately
    glafic_chi2_values = []  # Store glafic chi² for reference
    data_points_used = []  # Store number of data points used
    free_params_used = []  # Store number of free parameters used

    for output_file, fits_file, optresult_file in zip(output_files, fits_files, optresult_files):
        # Process output and fits files
        result = process_file_group(output_file, fits_file, input_file, n)
        results.append(result)
        
        # Calculate BIC using manual chi² but track both values
        bic, glafic_chi2, n_data_points, n_free_params = compute_bic_with_tracking(result['manual_chi2'], optresult_file)
        bics.append(bic)
        glafic_chi2_values.append(glafic_chi2)
        data_points_used.append(n_data_points)
        free_params_used.append(n_free_params)
        
        # Create LaTeX row with additional information
        model_name, constraint_name = determine_model_and_constraint(output_file)
        # Format individual position deltas and flux ratio differences for LaTeX
        pos_deltas_str = " & ".join([f"{delta:.3f}" for delta in result['individual_position_deltas']])
        flux_diffs_str = " & ".join([f"{diff:.3f}" for diff in result['flux_ratio_differences']])
        
        latex_row = f"{model_name} & {constraint_name} & {pos_deltas_str} & {result['delta_rms']:.1f} & {flux_diffs_str} & {result['predicted_flux_ratios_rms']:.2f} & {bic:.1f} \\\\"
        latex_rows.append(latex_row)
    
    # Combine results into a DataFrame
    df_results = pd.DataFrame(results)
    df_results['bic'] = bics  # Add BIC values to the DataFrame
    df_results['glafic_chi2'] = glafic_chi2_values  # Add glafic chi² for reference
    df_results['n_data_points'] = data_points_used  # Add number of data points used
    df_results['n_free_params'] = free_params_used  # Add number of free parameters used
    
    return df_results, latex_rows

# Determine model and constraint names from the filename
def determine_model_and_constraint(file_name):
    model_map = {'S': 'SIE', 'P': 'POW', 'N': 'NFW'}
    
    base_name = os.path.basename(file_name)
    model_key = base_name[3]
    constraint_key = base_name[4:6] if base_name[4:6] in ['PF'] else base_name[4]
    
    model_name = model_map.get(model_key, '')
    if 'R' in base_name:
        model_name += ' + Shear'
    if 'G' in base_name:
        model_name += ' + G2'
    
    constraint_name = 'QSO Pos + FR' if constraint_key == 'PF' else 'QSO Pos'
    
    return model_name, constraint_name

# Example usage with filenames
Input_file = "WFI2026/SPFC/Wobs_point.dat"

Output_files = [
    "WFI2026/SPC/outSP_point.dat",
    "WFI2026/SPC/outSPR_point.dat",
    "WFI2026/SPFC/outSPF_point.dat",
    "WFI2026/SPFC/outSPFR_point.dat",
    "WFI2026/PPC/outPP_point.dat",
    "WFI2026/PPC/outPPR_point.dat",
    "WFI2026/PPFC/outPPF_point.dat",
    "WFI2026/PPFC/outPPFR_point.dat",
    "WFI2026/NPC/outNP_point.dat",
    "WFI2026/NPC/outNPR_point.dat",
    "WFI2026/NPFC/outNPF_point.dat",
    "WFI2026/NPFC/outNPFR_point.dat",
]

Fits_files = [
    "WFI2026/SPC/outSP_lens.fits",
    "WFI2026/SPC/outSPR_lens.fits",
    "WFI2026/SPFC/outSPF_lens.fits",
    "WFI2026/SPFC/outSPFR_lens.fits",
    "WFI2026/PPC/outPP_lens.fits",
    "WFI2026/PPC/outPPR_lens.fits",
    "WFI2026/PPFC/outPPF_lens.fits",
    "WFI2026/PPFC/outPPFR_lens.fits",
    "WFI2026/NPC/outNP_lens.fits",
    "WFI2026/NPC/outNPR_lens.fits",
    "WFI2026/NPFC/outNPF_lens.fits",
    "WFI2026/NPFC/outNPFR_lens.fits",
]

Optresult_files = [
    "WFI2026/SPC/outSP_optresult.dat",
    "WFI2026/SPC/outSPR_optresult.dat",
    "WFI2026/SPFC/outSPF_optresult.dat",
    "WFI2026/SPFC/outSPFR_optresult.dat",
    "WFI2026/PPC/outPP_optresult.dat",
    "WFI2026/PPC/outPPR_optresult.dat",
    "WFI2026/PPFC/outPPF_optresult.dat",
    "WFI2026/PPFC/outPPFR_optresult.dat",
    "WFI2026/NPC/outNP_optresult.dat",
    "WFI2026/NPC/outNPR_optresult.dat",
    "WFI2026/NPFC/outNPF_optresult.dat",
    "WFI2026/NPFC/outNPFR_optresult.dat",
]

n = 0.0025  # Example scale factor
results_df, latex_rows = process_files(Input_file, Output_files, Fits_files, Optresult_files, n)

# Print formatted results
for i, output_file in enumerate(Output_files):
    constraint_type = "PF" if "PF" in Optresult_files[i] else "P"
    print(f"Model: {output_file}")
    print(f"Constraint Type: {constraint_type}")
    print(f"Number of data points used: {results_df.iloc[i]['n_data_points']}")
    print(f"Number of free parameters: {results_df.iloc[i]['n_free_params']}")
    print(f"Delta RMS of Image Positions: {results_df.iloc[i]['delta_rms']:.3f}")
    
    # Print individual position deltas (4 positions)
    print(f"Individual Position Deltas (mas):")
    for j, delta in enumerate(results_df.iloc[i]['individual_position_deltas']):
        print(f"  Image {j+1}: {delta:.3f}")
    
    print(f"Delta RMS of flux ratios at the observed position: {results_df.iloc[i]['predicted_flux_ratios_observed_rms']:.3f}")
    print(f"Delta RMS of flux ratios at the predicted position: {results_df.iloc[i]['predicted_flux_ratios_rms']:.3f}")
    
    # Print absolute flux ratio differences (3 ratios, excluding normalization)
    print(f"Absolute Flux Ratio Differences:")
    flux_diff_indices = [0, 2, 3] if len(results_df.iloc[i]['flux_ratio_differences']) == 3 else range(len(results_df.iloc[i]['flux_ratio_differences']))
    for j, diff in enumerate(results_df.iloc[i]['flux_ratio_differences']):
        img_idx = flux_diff_indices[j] + 1
        print(f"  Image {img_idx}: {diff:.3f}")
    
    print(f"Manual Chi²: {results_df.iloc[i]['manual_chi2']:.3f}")
    if results_df.iloc[i]['glafic_chi2'] is not None:
        print(f"Glafic Chi² (final, reference): {results_df.iloc[i]['glafic_chi2']:.3f}")
        print(f"Chi² Difference (Manual - Glafic): {results_df.iloc[i]['manual_chi2'] - results_df.iloc[i]['glafic_chi2']:.3f}")
    else:
        print(f"Glafic Chi² (reference): Not found")
    print(f"BIC (using manual Chi²): {results_df.iloc[i]['bic']:.3f}")
    
    # Manual verification of BIC calculation
    manual_bic_check = results_df.iloc[i]['manual_chi2'] + results_df.iloc[i]['n_free_params'] * np.log(results_df.iloc[i]['n_data_points'])
    print(f"Manual BIC verification: {results_df.iloc[i]['manual_chi2']:.3f} + {results_df.iloc[i]['n_free_params']} × ln({results_df.iloc[i]['n_data_points']}) = {manual_bic_check:.3f}")
    print()

# Print LaTeX-formatted rows
print("LaTeX Table Rows:")
print("Format: Model & Constraint & Pos1 & Pos2 & Pos3 & Pos4 & Overall ΔRMS & FluxDiff1 & FluxDiff2 & FluxDiff3 & Flux RMS & BIC")
for row in latex_rows:
    print(row)


# In[8]:


# NEW LATEX AUTOMATED CODE WITH MANUAL CHI^2 CALCULATION with individual position delta RMS and flux ratio differences
# WGDJ0405
import numpy as np
import pandas as pd
from astropy.io import fits
import os

# Function to calculate distance
def calculate_distance(x1, y1, x2, y2):
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

# Function to read observed positions, flux ratios, and uncertainties
def read_observed_positions(input_file):
    observed_positions = []
    observed_flux_ratios = []
    position_uncertainties = []
    flux_uncertainties = []
    with open(input_file, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line
        for line in lines:
            columns = line.split()
            if len(columns) >= 5:
                x1 = float(columns[0])
                y1 = float(columns[1])
                flux_ratio = float(columns[2])
                pos_err = float(columns[3])
                flux_err = float(columns[4])
                observed_positions.append((x1, y1))
                observed_flux_ratios.append(flux_ratio)
                position_uncertainties.append(pos_err)
                flux_uncertainties.append(flux_err)
    return observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties

# Function to read predicted positions and magnifications
def read_predicted_positions(output_file):
    predicted_positions = []
    predicted_magnifications = []
    with open(output_file, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line
        for line in lines:
            if line.startswith("#"):
                continue
            columns = line.split()
            if len(columns) >= 3:
                x2 = float(columns[0])
                y2 = float(columns[1])
                M2 = abs(float(columns[2]))
                predicted_positions.append((x2, y2))
                predicted_magnifications.append(M2)
    return predicted_positions, predicted_magnifications

# Function to pair predicted and observed positions
def pair_positions(observed_positions, predicted_positions, predicted_magnifications):
    paired_positions = []
    for observed_position in observed_positions:
        min_distance = float('inf')
        closest_index = None
        for i, predicted_position in enumerate(predicted_positions):
            distance = calculate_distance(*observed_position, *predicted_position)
            if distance < min_distance:
                min_distance = distance
                closest_index = i
        paired_positions.append((observed_position, predicted_positions[closest_index], predicted_magnifications[closest_index]))
    return paired_positions

# Function to calculate manual chi² with proper uncertainties
def calculate_manual_chi2(observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties, 
                         predicted_positions, predicted_magnifications, green_index=1):
    # Pair positions
    paired_positions = pair_positions(observed_positions, predicted_positions, predicted_magnifications)
    
    # Order predicted positions and magnifications based on observed order
    predicted_positions_ordered = [pair[1] for pair in paired_positions]
    predicted_magnifications_ordered = [pair[2] for pair in paired_positions]
    
    # Calculate predicted flux ratios
    predicted_flux_ratios = [m / predicted_magnifications_ordered[green_index] for m in predicted_magnifications_ordered]
    
    chi2 = 0.0
    
    # Position contributions to chi²
    for i, (observed, predicted, _) in enumerate(paired_positions):
        dx = observed[0] - predicted[0]
        dy = observed[1] - predicted[1]
        pos_err = position_uncertainties[i]
        chi2 += (dx/pos_err)**2 + (dy/pos_err)**2
    
    # Flux ratio contributions to chi² (excluding the normalization image)
    for i, (pred_flux, obs_flux) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)):
        if i != green_index:  # Exclude the green image used for normalization
            flux_err = flux_uncertainties[i]
            chi2 += ((pred_flux - obs_flux)/flux_err)**2
    
    return chi2

# Function to extract the final chi^2 from optresult.dat file
def extract_glafic_chi2(optresult_path):
    chi2_value = None
    with open(optresult_path, 'r') as file:
        lines = file.readlines()
    
    # Find the last occurrence of chi^2 (final/best result)
    for line in lines:
        if line.startswith("chi^2"):
            chi2_value = float(line.split()[2])
    
    return chi2_value

# Function to process file group
def process_file_group(output_file, fits_file, input_file, n):
    observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties = read_observed_positions(input_file)
    predicted_positions, predicted_magnifications = read_predicted_positions(output_file)
    
    paired_positions = pair_positions(observed_positions, predicted_positions, predicted_magnifications)
    
    # Calculate Delta RMS of Image Positions correctly
    squared_distances = 0.0
    individual_distances = []  # Store individual distances for each image
    for observed, predicted, _ in paired_positions:
        distance = calculate_distance(*observed, *predicted)
        individual_distances.append(distance * 1000)  # Convert to mas
        squared_distances += distance ** 2  # Summing squared distances
    
    mean_squared_distances = squared_distances / len(paired_positions)
    delta_rms = np.sqrt(mean_squared_distances)  # Taking square root of the mean of squared distances
    delta_rms *= 1000  # Convert to mas
    
    # Order predicted positions and magnifications based on observed order
    predicted_positions_ordered = [pair[1] for pair in paired_positions]
    predicted_magnifications_ordered = [pair[2] for pair in paired_positions]
    
    # Determine the index of the "Green Image"
    green_index = 2  # CHOOSE WHICH IMAGE TO NORMALISE BY HERE !!!!!!!!!!!!!!!!!!!!!!!

    # Predicted Flux Ratios
    predicted_flux_ratios = [m / predicted_magnifications_ordered[green_index] for m in predicted_magnifications_ordered]

    # Calculate absolute differences between predicted and observed flux ratios (excluding normalization image)
    flux_ratio_differences = []
    for i, (pred_flux, obs_flux) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)):
        if i != green_index:  # Exclude the green image used for normalization
            flux_ratio_differences.append(abs(pred_flux - obs_flux))

    # Calculate manual chi²
    manual_chi2 = calculate_manual_chi2(observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties,
                                       predicted_positions, predicted_magnifications, green_index)

    with fits.open(fits_file) as hdul:
        mag = hdul[0].data[6]
    
    pixel_coords = [(int((x + mag.shape[1] // 2 * n) / n), int((y + mag.shape[0] // 2 * n) / n)) for x, y in observed_positions]
    values = [mag[y, x] for x, y in pixel_coords]
    
    abs_values = [1 / abs(value) for value in values]
    predicted_flux_ratios_observed = [abs_values[i] / abs_values[green_index] for i in range(len(abs_values))]

    # Calculate RMS excluding the green image
    squared_diffs_pred = [(pred - obs) ** 2 for i, (pred, obs) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)) if i != green_index]
    predicted_flux_ratios_rms = np.sqrt(sum(squared_diffs_pred) / len(squared_diffs_pred))

    squared_diffs_obs = [(pred - obs) ** 2 for i, (pred, obs) in enumerate(zip(predicted_flux_ratios_observed, observed_flux_ratios)) if i != green_index]
    predicted_flux_ratios_observed_rms = np.sqrt(sum(squared_diffs_obs) / len(squared_diffs_obs))

    return {
        'delta_rms': delta_rms,
        'individual_position_deltas': individual_distances,  # Delta RMS for each of the 4 positions
        'flux_ratio_differences': flux_ratio_differences,    # Absolute differences for 3 flux ratios
        'predicted_flux_ratios_rms': predicted_flux_ratios_rms,
        'predicted_flux_ratios_observed_rms': predicted_flux_ratios_observed_rms,
        'manual_chi2': manual_chi2
    }

# Function to extract chi^2 and lens models from optresult.dat file (following Code 1's approach exactly)
def extract_chi2_and_lens_profiles(optresult_path):
    chi2_value = None
    profiles_used = []
    sie_count = 0

    with open(optresult_path, 'r') as file:
        lines = file.readlines()

    found_chi2 = False
    for line in lines:
        if line.startswith("chi^2"):
            chi2_value = float(line.split()[2])
            found_chi2 = True
            profiles_used = []
            sie_count = 0
        elif found_chi2:
            if line.startswith("lens   sie"):
                profiles_used.append("SIE")
                sie_count += 1
            elif line.startswith("lens   pow"):
                profiles_used.append("POW")
            elif line.startswith("lens   anfw"):
                profiles_used.append("NFW")
            elif line.startswith("lens   pert"):
                profiles_used.append("Shear")

    if chi2_value is None:
        raise ValueError("No chi^2 value found in the file.")

    profiles_used = list(set(profiles_used))
    return chi2_value, profiles_used, sie_count

# Function to calculate free parameters based on lens profiles (identical to Code 1)
def calculate_free_parameters(profiles_used, sie_count):
    profile_params = {
        "SIE": 5,
        "POW": 6,
        "NFW": 6,
        "Shear": 5
    }
    total_params = sum(profile_params[profile] for profile in profiles_used)
    if sie_count == 2:
        total_params += 5
    return total_params

# Function to calculate BIC (Bayesian Information Criterion)
def calculate_bic(chi2, n_data_points, n_free_params):
    return chi2 + n_free_params * np.log(n_data_points)

# Function to compute BIC using manual chi² but following Code 1's approach for free parameters
def compute_bic_with_tracking(manual_chi2, optresult_path):
    # Use Code 1's approach: extract profiles from optresult.dat file
    glafic_chi2, profiles_used, sie_count = extract_chi2_and_lens_profiles(optresult_path)
    n_free_params = calculate_free_parameters(profiles_used, sie_count)
    
    # Properly adjust n_data_points based on constraint type
    if "PF" in optresult_path:
        n_data_points = 4 * 2 + 3  # positions + flux ratios
    else:
        n_data_points = 4 * 2  # positions only
    
    # Always use manual chi² for BIC calculation
    bic_value = calculate_bic(manual_chi2, n_data_points, n_free_params)
    
    return bic_value, glafic_chi2, n_data_points, n_free_params

# Main function to process multiple files
def process_files(input_file, output_files, fits_files, optresult_files, n):
    results = []
    latex_rows = []
    bics = []  # Store BIC values separately
    glafic_chi2_values = []  # Store glafic chi² for reference
    data_points_used = []  # Store number of data points used
    free_params_used = []  # Store number of free parameters used

    for output_file, fits_file, optresult_file in zip(output_files, fits_files, optresult_files):
        # Process output and fits files
        result = process_file_group(output_file, fits_file, input_file, n)
        results.append(result)
        
        # Calculate BIC using manual chi² but track both values
        bic, glafic_chi2, n_data_points, n_free_params = compute_bic_with_tracking(result['manual_chi2'], optresult_file)
        bics.append(bic)
        glafic_chi2_values.append(glafic_chi2)
        data_points_used.append(n_data_points)
        free_params_used.append(n_free_params)
        
        # Create LaTeX row with additional information
        model_name, constraint_name = determine_model_and_constraint(output_file)
        # Format individual position deltas and flux ratio differences for LaTeX
        pos_deltas_str = " & ".join([f"{delta:.3f}" for delta in result['individual_position_deltas']])
        flux_diffs_str = " & ".join([f"{diff:.3f}" for diff in result['flux_ratio_differences']])
        
        latex_row = f"{model_name} & {constraint_name} & {pos_deltas_str} & {result['delta_rms']:.1f} & {flux_diffs_str} & {result['predicted_flux_ratios_rms']:.2f} & {bic:.1f} \\\\"
        latex_rows.append(latex_row)
    
    # Combine results into a DataFrame
    df_results = pd.DataFrame(results)
    df_results['bic'] = bics  # Add BIC values to the DataFrame
    df_results['glafic_chi2'] = glafic_chi2_values  # Add glafic chi² for reference
    df_results['n_data_points'] = data_points_used  # Add number of data points used
    df_results['n_free_params'] = free_params_used  # Add number of free parameters used
    
    return df_results, latex_rows

# Determine model and constraint names from the filename
def determine_model_and_constraint(file_name):
    model_map = {'S': 'SIE', 'P': 'POW', 'N': 'NFW'}
    
    base_name = os.path.basename(file_name)
    model_key = base_name[3]
    constraint_key = base_name[4:6] if base_name[4:6] in ['PF'] else base_name[4]
    
    model_name = model_map.get(model_key, '')
    if 'R' in base_name:
        model_name += ' + Shear'
    if 'G' in base_name:
        model_name += ' + G2'
    
    constraint_name = 'QSO Pos + FR' if constraint_key == 'PF' else 'QSO Pos'
    
    return model_name, constraint_name

# Example usage with filenames
Input_file = "WGDJ0405/SPFC/Eobs_point.dat"

Output_files = [
    "WGDJ0405/SPC/outSP_point.dat",
    "WGDJ0405/SPC/outSPR_point.dat",
    "WGDJ0405/SPFC/outSPF_point.dat",
    "WGDJ0405/SPFC/outSPFR_point.dat",
    "WGDJ0405/PPC/outPP_point.dat",
    "WGDJ0405/PPC/outPPR_point.dat",
    "WGDJ0405/PPFC/outPPF_point.dat",
    "WGDJ0405/PPFC/outPPFR_point.dat",
    "WGDJ0405/NPC/outNP_point.dat",
    "WGDJ0405/NPC/outNPR_point.dat",
    "WGDJ0405/NPFC/outNPF_point.dat",
    "WGDJ0405/NPFC/outNPFR_point.dat"
]

Fits_files = [
    "WGDJ0405/SPC/outSP_lens.fits",
    "WGDJ0405/SPC/outSPR_lens.fits",
    "WGDJ0405/SPFC/outSPF_lens.fits",
    "WGDJ0405/SPFC/outSPFR_lens.fits",
    "WGDJ0405/PPC/outPP_lens.fits",
    "WGDJ0405/PPC/outPPR_lens.fits",
    "WGDJ0405/PPFC/outPPF_lens.fits",
    "WGDJ0405/PPFC/outPPFR_lens.fits",
    "WGDJ0405/NPC/outNP_lens.fits",
    "WGDJ0405/NPC/outNPR_lens.fits",
    "WGDJ0405/NPFC/outNPF_lens.fits",
    "WGDJ0405/NPFC/outNPFR_lens.fits"
]

Optresult_files = [
    "WGDJ0405/SPC/outSP_optresult.dat",
    "WGDJ0405/SPC/outSPR_optresult.dat",
    "WGDJ0405/SPFC/outSPF_optresult.dat",
    "WGDJ0405/SPFC/outSPFR_optresult.dat",
    "WGDJ0405/PPC/outPP_optresult.dat",
    "WGDJ0405/PPC/outPPR_optresult.dat",
    "WGDJ0405/PPFC/outPPF_optresult.dat",
    "WGDJ0405/PPFC/outPPFR_optresult.dat",
    "WGDJ0405/NPC/outNP_optresult.dat",
    "WGDJ0405/NPC/outNPR_optresult.dat",
    "WGDJ0405/NPFC/outNPF_optresult.dat",
    "WGDJ0405/NPFC/outNPFR_optresult.dat"
]

n = 0.0025  # Example scale factor
results_df, latex_rows = process_files(Input_file, Output_files, Fits_files, Optresult_files, n)

# Print formatted results
for i, output_file in enumerate(Output_files):
    constraint_type = "PF" if "PF" in Optresult_files[i] else "P"
    print(f"Model: {output_file}")
    print(f"Constraint Type: {constraint_type}")
    print(f"Number of data points used: {results_df.iloc[i]['n_data_points']}")
    print(f"Number of free parameters: {results_df.iloc[i]['n_free_params']}")
    print(f"Delta RMS of Image Positions: {results_df.iloc[i]['delta_rms']:.3f}")
    
    # Print individual position deltas (4 positions)
    print(f"Individual Position Deltas (mas):")
    for j, delta in enumerate(results_df.iloc[i]['individual_position_deltas']):
        print(f"  Image {j+1}: {delta:.3f}")
    
    print(f"Delta RMS of flux ratios at the observed position: {results_df.iloc[i]['predicted_flux_ratios_observed_rms']:.3f}")
    print(f"Delta RMS of flux ratios at the predicted position: {results_df.iloc[i]['predicted_flux_ratios_rms']:.3f}")
    
    # Print absolute flux ratio differences (3 ratios, excluding normalization)
    print(f"Absolute Flux Ratio Differences:")
    flux_diff_indices = [0, 2, 3] if len(results_df.iloc[i]['flux_ratio_differences']) == 3 else range(len(results_df.iloc[i]['flux_ratio_differences']))
    for j, diff in enumerate(results_df.iloc[i]['flux_ratio_differences']):
        img_idx = flux_diff_indices[j] + 1
        print(f"  Image {img_idx}: {diff:.3f}")
    
    print(f"Manual Chi²: {results_df.iloc[i]['manual_chi2']:.3f}")
    if results_df.iloc[i]['glafic_chi2'] is not None:
        print(f"Glafic Chi² (final, reference): {results_df.iloc[i]['glafic_chi2']:.3f}")
        print(f"Chi² Difference (Manual - Glafic): {results_df.iloc[i]['manual_chi2'] - results_df.iloc[i]['glafic_chi2']:.3f}")
    else:
        print(f"Glafic Chi² (reference): Not found")
    print(f"BIC (using manual Chi²): {results_df.iloc[i]['bic']:.3f}")
    
    # Manual verification of BIC calculation
    manual_bic_check = results_df.iloc[i]['manual_chi2'] + results_df.iloc[i]['n_free_params'] * np.log(results_df.iloc[i]['n_data_points'])
    print(f"Manual BIC verification: {results_df.iloc[i]['manual_chi2']:.3f} + {results_df.iloc[i]['n_free_params']} × ln({results_df.iloc[i]['n_data_points']}) = {manual_bic_check:.3f}")
    print()

# Print LaTeX-formatted rows
print("LaTeX Table Rows:")
print("Format: Model & Constraint & Pos1 & Pos2 & Pos3 & Pos4 & Overall ΔRMS & FluxDiff1 & FluxDiff2 & FluxDiff3 & Flux RMS & BIC")
for row in latex_rows:
    print(row)


# In[12]:


# NEW LATEX AUTOMATED CODE WITH MANUAL CHI^2 CALCULATION with individual position delta RMS and flux ratio differences
# WGD2038
import numpy as np
import pandas as pd
from astropy.io import fits
import os

# Function to calculate distance
def calculate_distance(x1, y1, x2, y2):
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

# Function to read observed positions, flux ratios, and uncertainties
def read_observed_positions(input_file):
    observed_positions = []
    observed_flux_ratios = []
    position_uncertainties = []
    flux_uncertainties = []
    with open(input_file, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line
        for line in lines:
            columns = line.split()
            if len(columns) >= 5:
                x1 = float(columns[0])
                y1 = float(columns[1])
                flux_ratio = float(columns[2])
                pos_err = float(columns[3])
                flux_err = float(columns[4])
                observed_positions.append((x1, y1))
                observed_flux_ratios.append(flux_ratio)
                position_uncertainties.append(pos_err)
                flux_uncertainties.append(flux_err)
    return observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties

# Function to read predicted positions and magnifications
def read_predicted_positions(output_file):
    predicted_positions = []
    predicted_magnifications = []
    with open(output_file, 'r') as file:
        lines = file.readlines()[1:]  # Skip the first line
        for line in lines:
            if line.startswith("#"):
                continue
            columns = line.split()
            if len(columns) >= 3:
                x2 = float(columns[0])
                y2 = float(columns[1])
                M2 = abs(float(columns[2]))
                predicted_positions.append((x2, y2))
                predicted_magnifications.append(M2)
    return predicted_positions, predicted_magnifications

# Function to pair predicted and observed positions
def pair_positions(observed_positions, predicted_positions, predicted_magnifications):
    paired_positions = []
    for observed_position in observed_positions:
        min_distance = float('inf')
        closest_index = None
        for i, predicted_position in enumerate(predicted_positions):
            distance = calculate_distance(*observed_position, *predicted_position)
            if distance < min_distance:
                min_distance = distance
                closest_index = i
        paired_positions.append((observed_position, predicted_positions[closest_index], predicted_magnifications[closest_index]))
    return paired_positions

# Function to calculate manual chi² with proper uncertainties
def calculate_manual_chi2(observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties, 
                         predicted_positions, predicted_magnifications, green_index=1):
    # Pair positions
    paired_positions = pair_positions(observed_positions, predicted_positions, predicted_magnifications)
    
    # Order predicted positions and magnifications based on observed order
    predicted_positions_ordered = [pair[1] for pair in paired_positions]
    predicted_magnifications_ordered = [pair[2] for pair in paired_positions]
    
    # Calculate predicted flux ratios
    predicted_flux_ratios = [m / predicted_magnifications_ordered[green_index] for m in predicted_magnifications_ordered]
    
    chi2 = 0.0
    
    # Position contributions to chi²
    for i, (observed, predicted, _) in enumerate(paired_positions):
        dx = observed[0] - predicted[0]
        dy = observed[1] - predicted[1]
        pos_err = position_uncertainties[i]
        chi2 += (dx/pos_err)**2 + (dy/pos_err)**2
    
    # Flux ratio contributions to chi² (excluding the normalization image)
    for i, (pred_flux, obs_flux) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)):
        if i != green_index:  # Exclude the green image used for normalization
            flux_err = flux_uncertainties[i]
            chi2 += ((pred_flux - obs_flux)/flux_err)**2
    
    return chi2

# Function to extract the final chi^2 from optresult.dat file
def extract_glafic_chi2(optresult_path):
    chi2_value = None
    with open(optresult_path, 'r') as file:
        lines = file.readlines()
    
    # Find the last occurrence of chi^2 (final/best result)
    for line in lines:
        if line.startswith("chi^2"):
            chi2_value = float(line.split()[2])
    
    return chi2_value

# Function to process file group
def process_file_group(output_file, fits_file, input_file, n):
    observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties = read_observed_positions(input_file)
    predicted_positions, predicted_magnifications = read_predicted_positions(output_file)
    
    paired_positions = pair_positions(observed_positions, predicted_positions, predicted_magnifications)
    
    # Calculate Delta RMS of Image Positions correctly
    squared_distances = 0.0
    individual_distances = []  # Store individual distances for each image
    for observed, predicted, _ in paired_positions:
        distance = calculate_distance(*observed, *predicted)
        individual_distances.append(distance * 1000)  # Convert to mas
        squared_distances += distance ** 2  # Summing squared distances
    
    mean_squared_distances = squared_distances / len(paired_positions)
    delta_rms = np.sqrt(mean_squared_distances)  # Taking square root of the mean of squared distances
    delta_rms *= 1000  # Convert to mas
    
    # Order predicted positions and magnifications based on observed order
    predicted_positions_ordered = [pair[1] for pair in paired_positions]
    predicted_magnifications_ordered = [pair[2] for pair in paired_positions]
    
    # Determine the index of the "Green Image"
    green_index = 1  # CHOOSE WHICH IMAGE TO NORMALISE BY HERE !!!!!!!!!!!!!!!!!!!!!!!

    # Predicted Flux Ratios
    predicted_flux_ratios = [m / predicted_magnifications_ordered[green_index] for m in predicted_magnifications_ordered]

    # Calculate absolute differences between predicted and observed flux ratios (excluding normalization image)
    flux_ratio_differences = []
    for i, (pred_flux, obs_flux) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)):
        if i != green_index:  # Exclude the green image used for normalization
            flux_ratio_differences.append(abs(pred_flux - obs_flux))

    # Calculate manual chi²
    manual_chi2 = calculate_manual_chi2(observed_positions, observed_flux_ratios, position_uncertainties, flux_uncertainties,
                                       predicted_positions, predicted_magnifications, green_index)

    with fits.open(fits_file) as hdul:
        mag = hdul[0].data[6]
    
    pixel_coords = [(int((x + mag.shape[1] // 2 * n) / n), int((y + mag.shape[0] // 2 * n) / n)) for x, y in observed_positions]
    values = [mag[y, x] for x, y in pixel_coords]
    
    abs_values = [1 / abs(value) for value in values]
    predicted_flux_ratios_observed = [abs_values[i] / abs_values[green_index] for i in range(len(abs_values))]

    # Calculate RMS excluding the green image
    squared_diffs_pred = [(pred - obs) ** 2 for i, (pred, obs) in enumerate(zip(predicted_flux_ratios, observed_flux_ratios)) if i != green_index]
    predicted_flux_ratios_rms = np.sqrt(sum(squared_diffs_pred) / len(squared_diffs_pred))

    squared_diffs_obs = [(pred - obs) ** 2 for i, (pred, obs) in enumerate(zip(predicted_flux_ratios_observed, observed_flux_ratios)) if i != green_index]
    predicted_flux_ratios_observed_rms = np.sqrt(sum(squared_diffs_obs) / len(squared_diffs_obs))

    return {
        'delta_rms': delta_rms,
        'individual_position_deltas': individual_distances,  # Delta RMS for each of the 4 positions
        'flux_ratio_differences': flux_ratio_differences,    # Absolute differences for 3 flux ratios
        'predicted_flux_ratios_rms': predicted_flux_ratios_rms,
        'predicted_flux_ratios_observed_rms': predicted_flux_ratios_observed_rms,
        'manual_chi2': manual_chi2
    }

# Function to extract chi^2 and lens models from optresult.dat file (following Code 1's approach exactly)
def extract_chi2_and_lens_profiles(optresult_path):
    chi2_value = None
    profiles_used = []
    sie_count = 0

    with open(optresult_path, 'r') as file:
        lines = file.readlines()

    found_chi2 = False
    for line in lines:
        if line.startswith("chi^2"):
            chi2_value = float(line.split()[2])
            found_chi2 = True
            profiles_used = []
            sie_count = 0
        elif found_chi2:
            if line.startswith("lens   sie"):
                profiles_used.append("SIE")
                sie_count += 1
            elif line.startswith("lens   pow"):
                profiles_used.append("POW")
            elif line.startswith("lens   anfw"):
                profiles_used.append("NFW")
            elif line.startswith("lens   pert"):
                profiles_used.append("Shear")

    if chi2_value is None:
        raise ValueError("No chi^2 value found in the file.")

    profiles_used = list(set(profiles_used))
    return chi2_value, profiles_used, sie_count

# Function to calculate free parameters based on lens profiles (identical to Code 1)
def calculate_free_parameters(profiles_used, sie_count):
    profile_params = {
        "SIE": 5,
        "POW": 6,
        "NFW": 6,
        "Shear": 5
    }
    total_params = sum(profile_params[profile] for profile in profiles_used)
    if sie_count == 2:
        total_params += 5
    return total_params

# Function to calculate BIC (Bayesian Information Criterion)
def calculate_bic(chi2, n_data_points, n_free_params):
    return chi2 + n_free_params * np.log(n_data_points)

# Function to compute BIC using manual chi² but following Code 1's approach for free parameters
def compute_bic_with_tracking(manual_chi2, optresult_path):
    # Use Code 1's approach: extract profiles from optresult.dat file
    glafic_chi2, profiles_used, sie_count = extract_chi2_and_lens_profiles(optresult_path)
    n_free_params = calculate_free_parameters(profiles_used, sie_count)
    
    # Properly adjust n_data_points based on constraint type
    if "PF" in optresult_path:
        n_data_points = 4 * 2 + 3  # positions + flux ratios
    else:
        n_data_points = 4 * 2  # positions only
    
    # Always use manual chi² for BIC calculation
    bic_value = calculate_bic(manual_chi2, n_data_points, n_free_params)
    
    return bic_value, glafic_chi2, n_data_points, n_free_params

# Main function to process multiple files
def process_files(input_file, output_files, fits_files, optresult_files, n):
    results = []
    latex_rows = []
    bics = []  # Store BIC values separately
    glafic_chi2_values = []  # Store glafic chi² for reference
    data_points_used = []  # Store number of data points used
    free_params_used = []  # Store number of free parameters used

    for output_file, fits_file, optresult_file in zip(output_files, fits_files, optresult_files):
        # Process output and fits files
        result = process_file_group(output_file, fits_file, input_file, n)
        results.append(result)
        
        # Calculate BIC using manual chi² but track both values
        bic, glafic_chi2, n_data_points, n_free_params = compute_bic_with_tracking(result['manual_chi2'], optresult_file)
        bics.append(bic)
        glafic_chi2_values.append(glafic_chi2)
        data_points_used.append(n_data_points)
        free_params_used.append(n_free_params)
        
        # Create LaTeX row with additional information
        model_name, constraint_name = determine_model_and_constraint(output_file)
        # Format individual position deltas and flux ratio differences for LaTeX
        pos_deltas_str = " & ".join([f"{delta:.3f}" for delta in result['individual_position_deltas']])
        flux_diffs_str = " & ".join([f"{diff:.3f}" for diff in result['flux_ratio_differences']])
        
        latex_row = f"{model_name} & {constraint_name} & {pos_deltas_str} & {result['delta_rms']:.1f} & {flux_diffs_str} & {result['predicted_flux_ratios_rms']:.2f} & {bic:.1f} \\\\"
        latex_rows.append(latex_row)
    
    # Combine results into a DataFrame
    df_results = pd.DataFrame(results)
    df_results['bic'] = bics  # Add BIC values to the DataFrame
    df_results['glafic_chi2'] = glafic_chi2_values  # Add glafic chi² for reference
    df_results['n_data_points'] = data_points_used  # Add number of data points used
    df_results['n_free_params'] = free_params_used  # Add number of free parameters used
    
    return df_results, latex_rows

# Determine model and constraint names from the filename
def determine_model_and_constraint(file_name):
    model_map = {'S': 'SIE', 'P': 'POW', 'N': 'NFW'}
    
    base_name = os.path.basename(file_name)
    model_key = base_name[3]
    constraint_key = base_name[4:6] if base_name[4:6] in ['PF'] else base_name[4]
    
    model_name = model_map.get(model_key, '')
    if 'R' in base_name:
        model_name += ' + Shear'
    if 'G' in base_name:
        model_name += ' + G2'
    
    constraint_name = 'QSO Pos + FR' if constraint_key == 'PF' else 'QSO Pos'
    
    return model_name, constraint_name

# Example usage with filenames
Input_file = "WGD2038/SPFC/Robs_point_MT.dat"

Output_files = [
#    "WGD2038/SPC/outSP_point.dat",
#    "WGD2038/SPC/outSPR_point.dat",
#    "WGD2038/SPFC/outSPF_point.dat",
#    "WGD2038/SPFC/outSPFR_point.dat",
    "WGD2038/SPFC/outSPF_MT_point.dat",
#    "WGD2038/SPFC/outSPF_PT_point.dat",
#    "WGD2038/PPC/outPP_point.dat",
#    "WGD2038/PPC/outPPR_point.dat",
#    "WGD2038/PPFC/outPPF_point.dat",
#    "WGD2038/PPFC/outPPFR_point.dat",
#    "WGD2038/NPC/outNP_point.dat",
#    "WGD2038/NPC/outNPR_point.dat",
#    "WGD2038/NPFC/outNPF_point.dat",
#    "WGD2038/NPFC/outNPFR_point.dat"
]

Fits_files = [
#    "WGD2038/SPC/outSP_lens.fits",
#    "WGD2038/SPC/outSPR_lens.fits",
#    "WGD2038/SPFC/outSPF_lens.fits",
#    "WGD2038/SPFC/outSPFR_lens.fits",
    "WGD2038/SPFC/outSPF_MT_lens.fits",
#    "WGD2038/SPFC/outSPF_PT_lens.fits",
#    "WGD2038/PPC/outPP_lens.fits",
#    "WGD2038/PPC/outPPR_lens.fits",
#    "WGD2038/PPFC/outPPF_lens.fits",
#    "WGD2038/PPFC/outPPFR_lens.fits",
#    "WGD2038/NPC/outNP_lens.fits",
#    "WGD2038/NPC/outNPR_lens.fits",
#    "WGD2038/NPFC/outNPF_lens.fits",
#    "WGD2038/NPFC/outNPFR_lens.fits"
]

Optresult_files = [
#    "WGD2038/SPC/outSP_optresult.dat",
#    "WGD2038/SPC/outSPR_optresult.dat",
#    "WGD2038/SPFC/outSPF_optresult.dat",
#    "WGD2038/SPFC/outSPFR_optresult.dat",
    "WGD2038/SPFC/outSPF_MT_optresult.dat",
#    "WGD2038/SPFC/outSPF_PT_optresult.dat",
#    "WGD2038/PPC/outPP_optresult.dat",
#    "WGD2038/PPC/outPPR_optresult.dat",
#    "WGD2038/PPFC/outPPF_optresult.dat",
#    "WGD2038/PPFC/outPPFR_optresult.dat",
#    "WGD2038/NPC/outNP_optresult.dat",
#    "WGD2038/NPC/outNPR_optresult.dat",
#    "WGD2038/NPFC/outNPF_optresult.dat",
#    "WGD2038/NPFC/outNPFR_optresult.dat"
]

n = 0.0025  # Example scale factor
results_df, latex_rows = process_files(Input_file, Output_files, Fits_files, Optresult_files, n)

# Print formatted results
for i, output_file in enumerate(Output_files):
    constraint_type = "PF" if "PF" in Optresult_files[i] else "P"
    print(f"Model: {output_file}")
    print(f"Constraint Type: {constraint_type}")
    print(f"Number of data points used: {results_df.iloc[i]['n_data_points']}")
    print(f"Number of free parameters: {results_df.iloc[i]['n_free_params']}")
    print(f"Delta RMS of Image Positions: {results_df.iloc[i]['delta_rms']:.3f}")
    
    # Print individual position deltas (4 positions)
    print(f"Individual Position Deltas (mas):")
    for j, delta in enumerate(results_df.iloc[i]['individual_position_deltas']):
        print(f"  Image {j+1}: {delta:.3f}")
    
    print(f"Delta RMS of flux ratios at the observed position: {results_df.iloc[i]['predicted_flux_ratios_observed_rms']:.3f}")
    print(f"Delta RMS of flux ratios at the predicted position: {results_df.iloc[i]['predicted_flux_ratios_rms']:.3f}")
    
    # Print absolute flux ratio differences (3 ratios, excluding normalization)
    print(f"Absolute Flux Ratio Differences:")
    flux_diff_indices = [0, 2, 3] if len(results_df.iloc[i]['flux_ratio_differences']) == 3 else range(len(results_df.iloc[i]['flux_ratio_differences']))
    for j, diff in enumerate(results_df.iloc[i]['flux_ratio_differences']):
        img_idx = flux_diff_indices[j] + 1
        print(f"  Image {img_idx}: {diff:.3f}")
    
    print(f"Manual Chi²: {results_df.iloc[i]['manual_chi2']:.3f}")
    if results_df.iloc[i]['glafic_chi2'] is not None:
        print(f"Glafic Chi² (final, reference): {results_df.iloc[i]['glafic_chi2']:.3f}")
        print(f"Chi² Difference (Manual - Glafic): {results_df.iloc[i]['manual_chi2'] - results_df.iloc[i]['glafic_chi2']:.3f}")
    else:
        print(f"Glafic Chi² (reference): Not found")
    print(f"BIC (using manual Chi²): {results_df.iloc[i]['bic']:.3f}")
    
    # Manual verification of BIC calculation
    manual_bic_check = results_df.iloc[i]['manual_chi2'] + results_df.iloc[i]['n_free_params'] * np.log(results_df.iloc[i]['n_data_points'])
    print(f"Manual BIC verification: {results_df.iloc[i]['manual_chi2']:.3f} + {results_df.iloc[i]['n_free_params']} × ln({results_df.iloc[i]['n_data_points']}) = {manual_bic_check:.3f}")
    print()

# Print LaTeX-formatted rows
print("LaTeX Table Rows:")
print("Format: Model & Constraint & Pos1 & Pos2 & Pos3 & Pos4 & Overall ΔRMS & FluxDiff1 & FluxDiff2 & FluxDiff3 & Flux RMS & BIC")
for row in latex_rows:
    print(row)

