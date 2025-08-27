#!/usr/bin/env python
# coding: utf-8

# In[2]:


# Python script to process and rename optresult files and save the output to a text file


def extract_shear_and_convergence(filename):
    """
    Extract the shear amplitude and convergence values from the latest "lens   pert" line in the file.
    """
    
    shear_amplitude = None
    convergence = None

    with open(filename, 'r') as file:
        lines = file.readlines()

    # Iterate through the file from bottom to top to find the latest "lens   pert" line
    for line in reversed(lines):
        if "lens   sie" in line:
            # Extract the numerical values correctly using a fixed format approach
            # Split the line by whitespace but ensure consistent indexing
            parts = [value for value in line.split() if value]

            try:
                # Extract the 5th and 8th numbers (counted left-to-right)
                shear_amplitude = float(parts[6])  # 5th parameter
                convergence = float(parts[9])     # 8th parameter

            except (ValueError, IndexError):
                print(f"Error: Failed to extract numbers from file: {filename}. Check the 'lens pert' line format.")
                return None, None

            # Once found, stop searching
            break

    if shear_amplitude is not None and convergence is not None:
        return shear_amplitude, convergence

    else:
        print(f"Error: 'lens pert' line not found in file: {filename}.")
        return None, None


def rename_file(file_path):
    """
    Rename a file based on its naming conventions.
    """
    
    # Split the file path to extract components
    parts = file_path.split("/")
    lens_system = parts[0]  # Extract the lens system (e.g., RXJ0911)
    file_name = parts[-1]  # Extract the file name (e.g., outSPFR_optresult.dat)

    # Extract the relevant part of the filename before "_optresult"
    key = file_name.split("_optresult")[0][3:]  # Remove "out" prefix

    # Determine the profile type from the first letter
    profile_map = {"S": "SIE", "P": "POW", "N": "NFW"}
    profile = profile_map.get(key[0], "Unknown")

    # Determine the constraint type based on the second character or characters
    if key[1:3] == "PF":
        constraint = "position and flux constraint"

    elif key[1] == "P":
        constraint = "position constraint"

    else:
        constraint = "Unknown constraint"

    # Determine if "G2" or "Shear" is present
    extras = []
    if "G" in key:
        extras.append("G2")

    if "R" in key:
        extras.append("Shear")

    # Construct the new name
    new_name = f"{lens_system}: {profile} {constraint}"
    if extras:
        new_name += ", " + " + ".join(extras)

    return new_name


def process_multiple_files(file_list, output_file):
    """
    Process a list of files, rename them, and write the extracted results to an output file.
    """
    
    with open(output_file, 'w') as out_file:
        for file_name in file_list:
            # Rename the file based on naming conventions
            description = rename_file(file_name)

            # Extract shear amplitude and convergence
            shear_amplitude, convergence = extract_shear_and_convergence(file_name)

            if shear_amplitude is not None and convergence is not None:
                # Format the output and write to the file
                output_line = f"{description} : Shear Amplitude = {shear_amplitude}, Convergence = {convergence}\n\n"
                out_file.write(output_line)
                print(output_line.strip())  # Print to console without the trailing newline

            else:
                # Log an error entry if values could not be extracted
                error_line = f"{description}: Failed to extract shear amplitude and convergence.\n\n"
                out_file.write(error_line)
                print(error_line.strip())  # Print to console without the trailing newline


# List of file paths to process
file_list = [
    "RXJ0911/SPC/outSP_optresult.dat",
    "RXJ0911/SPC/outSPR_optresult.dat",
    "RXJ0911/SPC/outSPG_optresult.dat",
    "RXJ0911/SPC/outSPGR_optresult.dat",
    "RXJ0911/SPFC/outSPF_optresult.dat",
    "RXJ0911/SPFC/outSPFR_optresult.dat",
    "RXJ0911/SPFC/outSPFG_optresult.dat",
    "RXJ0911/SPFC/outSPFGR_optresult.dat",
]

# Output file to store results
output_file = "SIE_Parameter_results_test.txt"

# Process the files and save output to the specified file
process_multiple_files(file_list, output_file)


# In[4]:


def extract_sie_parameters(filename):
    """
    Extract the SIE lens parameters from the latest "lens   sie" line in the file.
    """
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Iterate through the file from bottom to top to find the latest "lens   sie" line
    for line in reversed(lines):
        if "lens   sie" in line:
            # Split the line by whitespace but ensure consistent indexing
            parts = line.split()

            try:
                # Extract the specified parameters using index position
                lens_redshift = float(parts[2])
                velocity_dispersion = float(parts[3])
                x_coord = float(parts[4])
                y_coord = float(parts[5])
                ellipticity = float(parts[6])
                position_angle = float(parts[7])

            except (ValueError, IndexError):
                print(f"Error: Failed to extract numbers from file: {filename}. Check the 'lens sie' line format.")
                return None

            # Once found, return the extracted values as a tuple
            return (lens_redshift, velocity_dispersion, x_coord, y_coord, ellipticity, position_angle)

    print(f"Error: 'lens sie' line not found in file: {filename}.")
    return None


def process_multiple_files(file_list, output_file):
    """
    Process a list of files, rename them, and write the extracted SIE lens parameters to an output file.
    """
    with open(output_file, 'w') as out_file:
        for file_name in file_list:
            # Rename the file based on naming conventions
            description = rename_file(file_name)

            # Extract SIE lens parameters
            sie_parameters = extract_sie_parameters(file_name)

            if sie_parameters:
                # Format the output and write to the file
                output_line = f"{description} : Lens Redshift = {sie_parameters[0]}, Velocity Dispersion = {sie_parameters[1]}, X-coord = {sie_parameters[2]}, Y-coord = {sie_parameters[3]}, Ellipticity = {sie_parameters[4]}, Position Angle = {sie_parameters[5]}\n\n"
                out_file.write(output_line)
                print(output_line.strip())  # Print to console without the trailing newline

            else:
                # Log an error entry if values could not be extracted
                error_line = f"{description}: Failed to extract SIE parameters.\n\n"
                out_file.write(error_line)
                print(error_line.strip())  # Print to console without the trailing newline


# List of file paths to process
file_list = [
    "RXJ0911/SPC/outSP_optresult.dat",
    "RXJ0911/SPC/outSPR_optresult.dat",
    "RXJ0911/SPC/outSPG_optresult.dat",
    "RXJ0911/SPC/outSPGR_optresult.dat",
    "RXJ0911/SPFC/outSPF_optresult.dat",
    "RXJ0911/SPFC/outSPFR_optresult.dat",
    "RXJ0911/SPFC/outSPFG_optresult.dat",
    "RXJ0911/SPFC/outSPFGR_optresult.dat",
]

# Output file to store results
output_file = "SIE_Parameter_results_test.txt"

# Process the files and save output to the specified file
process_multiple_files(file_list, output_file)


# In[8]:


def extract_sie_parameters(filename):
    """
    Extract the SIE lens parameters from the first "lens   sie" line that appears after the final mention of "chi^2" in the file.
    """
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Locate the last occurrence of "chi^2"
    last_chi2_index = None
    for index, line in enumerate(lines):
        if "chi^2" in line:
            last_chi2_index = index

    # If no "chi^2" was found, we can't process further
    if last_chi2_index is None:
        print(f"Error: 'chi^2' line not found in file: {filename}.")
        return None

    # Scan downwards from the last "chi^2" line to find the first "lens sie" line
    for line in lines[last_chi2_index:]:
        if "lens   sie" in line:
            parts = line.split()
            try:
                # The index positions may need to be adjusted based on the specific formatting of the data
                lens_redshift = float(parts[2])
                velocity_dispersion = float(parts[3])
                x_coord = float(parts[4])
                y_coord = float(parts[5])
                ellipticity = float(parts[6])
                position_angle = float(parts[7])
                return (lens_redshift, velocity_dispersion, x_coord, y_coord, ellipticity, position_angle)
            except (ValueError, IndexError):
                print(f"Error: Failed to extract numbers from file: {filename}. Check the 'lens sie' line format.")
                return None

    print(f"Error: 'lens sie' line not found after last 'chi^2' in file: {filename}.")
    return None

def rename_file(file_path):
    """
    Rename a file based on its naming conventions.
    """
    parts = file_path.split("/")
    lens_system = parts[0]  # Extract the lens system
    file_name = parts[-1]  # Extract the file name

    key = file_name.split("_optresult")[0][3:]  # Remove "out" prefix
    profile_map = {"S": "SIE", "P": "POW", "N": "NFW"}
    profile = profile_map.get(key[0], "Unknown")

    constraint = "Unknown constraint"
    if key[1:3] == "PF":
        constraint = "position and flux constraint"
    elif key[1] == "P":
        constraint = "position constraint"

    extras = []
    if "G" in key:
        extras.append("G2")
    if "R" in key:
        extras.append("Shear")

    new_name = f"{lens_system}: {profile} {constraint}" + (", " + " + ".join(extras) if extras else "")
    return new_name

def process_multiple_files(file_list, output_file):
    """
    Process a list of files, rename them, and write the extracted SIE lens parameters to an output file.
    """
    with open(output_file, 'w') as out_file:
        for file_name in file_list:
            description = rename_file(file_name)
            sie_parameters = extract_sie_parameters(file_name)

            if sie_parameters:
                output_line = f"{description} : Lens Redshift = {sie_parameters[0]}, Velocity Dispersion = {sie_parameters[1]}, X-coord = {sie_parameters[2]}, Y-coord = {sie_parameters[3]}, Ellipticity = {sie_parameters[4]}, Position Angle = {sie_parameters[5]}\n\n"
                out_file.write(output_line)
                print(output_line.strip())
            else:
                error_line = f"{description}: Failed to extract SIE parameters.\n\n"
                out_file.write(error_line)
                print(error_line.strip())

# List of file paths to process
file_list = [
    "RXJ0911/SPC/outSP_optresult.dat",
    "RXJ0911/SPC/outSPR_optresult.dat",
    "RXJ0911/SPC/outSPG_optresult.dat",
    "RXJ0911/SPC/outSPGR_optresult.dat",
    "RXJ0911/SPFC/outSPF_optresult.dat",
    "RXJ0911/SPFC/outSPFR_optresult.dat",
    "RXJ0911/SPFC/outSPFG_optresult.dat",
    "RXJ0911/SPFC/outSPFGR_optresult.dat",
]

# Output file to store results
output_file = "SIE_Parameter_results_test.txt"

# Process the files and save output to the specified file
process_multiple_files(file_list, output_file)


# In[14]:


def extract_sie_parameters(filename):
    """
    Extract the SIE lens parameters from the first "lens   sie" line that appears after the final mention of "chi^2" in the file.
    """
    with open(filename, 'r') as file:
        lines = file.readlines()

    last_chi2_index = None
    for index, line in enumerate(lines):
        if "chi^2" in line:
            last_chi2_index = index

    if last_chi2_index is None:
        print(f"Error: 'chi^2' line not found in file: {filename}.")
        return None

    for line in lines[last_chi2_index:]:
        if "lens   sie" in line:
            parts = line.split()
            try:
                lens_redshift = float(parts[2])
                velocity_dispersion = float(parts[3])
                x_coord = float(parts[4])
                y_coord = float(parts[5])
                ellipticity = float(parts[6])
                position_angle = float(parts[7])
                return (lens_redshift, velocity_dispersion, x_coord, y_coord, ellipticity, position_angle)
            except (ValueError, IndexError):
                print(f"Error: Failed to extract numbers from file: {filename}. Check the 'lens sie' line format.")
                return None

    print(f"Error: 'lens sie' line not found after last 'chi^2' in file: {filename}.")
    return None

def rename_file(file_path):
    parts = file_path.split("/")
    lens_system = parts[0]
    file_name = parts[-1]

    key = file_name.split("_optresult")[0][3:]
    profile_map = {"S": "SIE", "P": "POW", "N": "NFW"}
    profile = profile_map.get(key[0], "Unknown")

    constraint = "Unknown constraint"
    if key[1:3] == "PF":
        constraint = "position and flux constraint"
    elif key[1] == "P":
        constraint = "position constraint"

    extras = []
    if "G" in key:
        extras.append("G2")
    if "R" in key:
        extras.append("Shear")

    new_name = f"{lens_system}: {profile} {constraint}" + (", " + " + ".join(extras) if extras else "")
    return new_name

def process_multiple_files(file_list):
    results = []
    for file_name in file_list:
        model_description = rename_file(file_name)
        sie_parameters = extract_sie_parameters(file_name)
        if sie_parameters:
            results.append(f"{model_description} & {sie_parameters[0]:.3f} & {sie_parameters[1]:.3f} & {sie_parameters[2]:.3f} & {sie_parameters[3]:.3f} & {sie_parameters[4]:.3f} & {sie_parameters[5]:.3f} \\\\")
        else:
            results.append(f"{model_description}: Failed to extract SIE parameters.")
    return results

# List of file paths to process
file_list = [
    "RXJ0911/SPC/outSP_optresult.dat",
    "RXJ0911/SPC/outSPR_optresult.dat",
    "RXJ0911/SPC/outSPG_optresult.dat",
    "RXJ0911/SPC/outSPGR_optresult.dat",
    "RXJ0911/SPFC/outSPF_optresult.dat",
    "RXJ0911/SPFC/outSPFR_optresult.dat",
    "RXJ0911/SPFC/outSPFG_optresult.dat",
    "RXJ0911/SPFC/outSPFGR_optresult.dat",
    "PSJ1606/SPC/outSP_optresult.dat"
]

# Process the files and print LaTeX-formatted results
latex_results = process_multiple_files(file_list)
for result in latex_results:
    print(result)


# In[21]:


def extract_sie_parameters(filename):
    """
    Extract the SIE lens parameters from the first "lens   sie" line that appears after the final mention of "chi^2" in the file.
    """
    with open(filename, 'r') as file:
        lines = file.readlines()

    last_chi2_index = None
    for index, line in enumerate(lines):
        if "chi^2" in line:
            last_chi2_index = index

    if last_chi2_index is None:
        print(f"Error: 'chi^2' line not found in file: {filename}.")
        return None

    for line in lines[last_chi2_index:]:
        if "lens   sie" in line:
            parts = line.split()
            try:
                lens_redshift = float(parts[2])
                velocity_dispersion = float(parts[3])
                x_coord = float(parts[4])
                y_coord = float(parts[5])
                ellipticity = float(parts[6])
                position_angle = float(parts[7])
                return (lens_redshift, velocity_dispersion, x_coord, y_coord, ellipticity, position_angle)
            except (ValueError, IndexError):
                print(f"Error: Failed to extract numbers from file: {filename}. Check the 'lens sie' line format.")
                return None

    print(f"Error: 'lens sie' line not found after last 'chi^2' in file: {filename}.")
    return None

def rename_file(file_path):
    parts = file_path.split("/")
    lens_system = parts[0]  # Now correctly fetching the first part of the file path
    file_name = parts[-1]

    key = file_name.split("_optresult")[0][3:]
    profile_map = {"S": "SIE", "P": "POW", "N": "NFW"}
    profile = profile_map.get(key[0], "Unknown")

    constraint = "Unknown constraint"
    if key[1:3] == "PF":
        constraint = "position and flux constraint"
    elif key[1] == "P":
        constraint = "position constraint"

    extras = []
    if "G" in key:
        extras.append("G2")
    if "R" in key:
        extras.append("Shear")

    new_name = f"{lens_system} & {profile} {constraint}" + (", " + " + ".join(extras) if extras else "")
    return new_name

def process_multiple_files(file_list):
    results = []
    for file_name in file_list:
        model_description = rename_file(file_name)
        sie_parameters = extract_sie_parameters(file_name)
        if sie_parameters:
            results.append(f"{model_description} & {sie_parameters[0]:.3f} & {sie_parameters[1]:.3f} & {sie_parameters[2]:.3f} & {sie_parameters[3]:.3f} & {sie_parameters[4]:.3f} & {sie_parameters[5]:.3f} \\\\")
        else:
            results.append(f"{model_description} Failed to extract SIE parameters.")
    return results

# List of file paths to process
file_list = [
    "RXJ0911/SPC/outSP_optresult.dat",
    "RXJ0911/SPC/outSPR_optresult.dat",
    "RXJ0911/SPC/outSPG_optresult.dat",
    "RXJ0911/SPC/outSPGR_optresult.dat",
    "RXJ0911/SPFC/outSPF_optresult.dat",
    "RXJ0911/SPFC/outSPFR_optresult.dat",
    "RXJ0911/SPFC/outSPFG_optresult.dat",
    "RXJ0911/SPFC/outSPFGR_optresult.dat",
    "PSJ1606/SPC/outSP_optresult.dat"
]

def process_multiple_files(file_list):
    results = []
    for file_name in file_list:
        model_description = rename_file(file_name)
        sie_parameters = extract_sie_parameters(file_name)
        if sie_parameters:
            lens_redshift, velocity_dispersion, x_coord, y_coord, ellipticity, position_angle = sie_parameters
            result_text = f"{model_description} & {lens_redshift:.3f} & {velocity_dispersion:.3f} & {x_coord:.3f} & {y_coord:.3f} & {ellipticity:.3f} & {position_angle:.3f} \\\\"
            results.append(result_text)
            print(f"Model: {file_name}")
            print(f"Lens Redshift: {lens_redshift:.3f}")
            print(f"Velocity Dispersion: {velocity_dispersion:.3f}")
            print(f"X Coordinate: {x_coord:.3f}")
            print(f"Y Coordinate: {y_coord:.3f}")
            print(f"Ellipticity: {ellipticity:.3f}")
            print(f"Position Angle: {position_angle:.3f}\n")
        else:
            results.append(f"{model_description} Failed to extract SIE parameters.")
            print(f"Model: {file_name} Failed to extract SIE parameters.\n")
    return results

# Process the files and print LaTeX-formatted results
latex_results = process_multiple_files(file_list)
for result in latex_results:
    print(result)


# In[23]:


def extract_sie_parameters(filename):
    """
    Extract the SIE lens parameters from the first "lens   sie" line that appears after the final mention of "chi^2" in the file.
    """
    with open(filename, 'r') as file:
        lines = file.readlines()

    last_chi2_index = None
    for index, line in enumerate(lines):
        if "chi^2" in line:
            last_chi2_index = index

    if last_chi2_index is None:
        print(f"Error: 'chi^2' line not found in file: {filename}.")
        return None

    for line in lines[last_chi2_index:]:
        if "lens   sie" in line:
            parts = line.split()
            try:
                lens_redshift = float(parts[2])
                velocity_dispersion = float(parts[3])
                x_coord = float(parts[4])
                y_coord = float(parts[5])
                ellipticity = float(parts[6])
                position_angle = float(parts[7])
                return (lens_redshift, velocity_dispersion, x_coord, y_coord, ellipticity, position_angle)
            except (ValueError, IndexError):
                print(f"Error: Failed to extract numbers from file: {filename}. Check the 'lens sie' line format.")
                return None

    print(f"Error: 'lens sie' line not found after last 'chi^2' in file: {filename}.")
    return None

def rename_file(file_path):
    parts = file_path.split("/")
    lens_system = parts[0]  # Now correctly fetching the first part of the file path
    file_name = parts[-1]

    key = file_name.split("_optresult")[0][3:]
    profile_map = {"S": "SIE", "P": "POW", "N": "NFW"}
    profile = profile_map.get(key[0], "Unknown")

    constraint = "Unknown constraint"
    if key[1:3] == "PF":
        constraint = "position and flux constraint"
    elif key[1] == "P":
        constraint = "position constraint"

    extras = []
    if "G" in key:
        extras.append("G2")
    if "R" in key:
        extras.append("Shear")

    new_name = f"{lens_system} & {profile} {constraint}" + (", " + " + ".join(extras) if extras else "")
    return new_name

def process_multiple_files(file_list):
    results = []
    for file_name in file_list:
        model_description = rename_file(file_name)
        sie_parameters = extract_sie_parameters(file_name)
        if sie_parameters:
            results.append(f"{model_description} & {sie_parameters[0]:.3f} & {sie_parameters[1]:.3f} & {sie_parameters[2]:.3f} & {sie_parameters[3]:.3f} & {sie_parameters[4]:.3f} & {sie_parameters[5]:.3f} \\\\")
        else:
            results.append(f"{model_description} Failed to extract SIE parameters.")
    return results

# List of file paths to process
file_list = [
    "RXJ0911/SPC/outSP_optresult.dat",
    "RXJ0911/SPC/outSPR_optresult.dat",
    "RXJ0911/SPC/outSPG_optresult.dat",
    "RXJ0911/SPC/outSPGR_optresult.dat",
    "RXJ0911/SPFC/outSPF_optresult.dat",
    "RXJ0911/SPFC/outSPFR_optresult.dat",
    "RXJ0911/SPFC/outSPFG_optresult.dat",
    "RXJ0911/SPFC/outSPFGR_optresult.dat",
    "PSJ1606/SPC/outSP_optresult.dat",
    "PSJ1606/SPC/outSPR_optresult.dat",
    "PSJ1606/SPC/outSPG_optresult.dat",
    "PSJ1606/SPC/outSPGR_optresult.dat",
    "PSJ1606/SPFC/outSPF_optresult.dat",
    "PSJ1606/SPFC/outSPFR_optresult.dat",
    "PSJ1606/SPFC/outSPFG_optresult.dat",
    "PSJ1606/SPFC/outSPFGR_optresult.dat",
    "WFI2033/SP/outSP_optresult.dat",
    "WFI2033/SPR/outSPR_optresult.dat",
    "WFI2033/SPG/outSPG_optresult.dat",
    "WFI2033/SPGR/outSPGR_optresult.dat",
    "WFI2033/SPF/outSPF_optresult.dat",
    "WFI2033/SPFR/outSPFR_optresult.dat",
    "WFI2033/SPFG/outSPFG_optresult.dat",
    "WFI2033/SPFGR/outSPFGR_optresult.dat",
    "SDSSJ1330/SPC/outSP_optresult.dat",
    "SDSSJ1330/SPC/outSPR_optresult.dat",
    "SDSSJ1330/SPFC/outSPF_optresult.dat",
    "SDSSJ1330/SPFC/outSPFR_optresult.dat",
    "WFI2026/SPC/outSP_optresult.dat",
    "WFI2026/SPC/outSPR_optresult.dat",
    "WFI2026/SPFC/outSPF_optresult.dat",
    "WFI2026/SPFC/outSPFR_optresult.dat",
    "WGDJ0405/SPC/outSP_optresult.dat",
    "WGDJ0405/SPC/outSPR_optresult.dat",
    "WGDJ0405/SPFC/outSPF_optresult.dat",
    "WGDJ0405/SPFC/outSPFR_optresult.dat"
]

def process_multiple_files(file_list):
    results = []
    for file_name in file_list:
        model_description = rename_file(file_name)
        sie_parameters = extract_sie_parameters(file_name)
        if sie_parameters:
            lens_redshift, velocity_dispersion, x_coord, y_coord, ellipticity, position_angle = sie_parameters
            result_text = f"{model_description} & {lens_redshift:.3f} & {velocity_dispersion:.3f} & {x_coord:.3f} & {y_coord:.3f} & {ellipticity:.3f} & {position_angle:.3f} \\\\"
            results.append(result_text)
            print(f"Model: {file_name}")
            print(f"Lens Redshift: {lens_redshift:.3f}")
            print(f"Velocity Dispersion: {velocity_dispersion:.3f}")
            print(f"X Coordinate: {x_coord:}")
            print(f"Y Coordinate: {y_coord:}")
            print(f"Ellipticity: {ellipticity:.3f}")
            print(f"Position Angle: {position_angle:.3f}\n")
        else:
            results.append(f"{model_description} Failed to extract SIE parameters.")
            print(f"Model: {file_name} Failed to extract SIE parameters.\n")
    return results

# Process the files and print LaTeX-formatted results
latex_results = process_multiple_files(file_list)
for result in latex_results:
    print(result)


# In[2]:


# Successful Code that outputs the parameters for SIE as required!!! 

def extract_sie_parameters(filename):
    """
    Extract the SIE lens parameters from the first "lens   sie" line that appears after the final mention of "chi^2" in the file.
    """
    with open(filename, 'r') as file:
        lines = file.readlines()

    last_chi2_index = None
    for index, line in enumerate(lines):
        if "chi^2" in line:
            last_chi2_index = index

    if last_chi2_index is None:
        print(f"Error: 'chi^2' line not found in file: {filename}.")
        return None

    for line in lines[last_chi2_index:]:
        if "lens   sie" in line:
            parts = line.split()
            try:
                lens_redshift = float(parts[2])
                velocity_dispersion = float(parts[3])
                x_coord = float(parts[4])
                y_coord = float(parts[5])
                ellipticity = float(parts[6])
                position_angle = float(parts[7])
                return (lens_redshift, velocity_dispersion, x_coord, y_coord, ellipticity, position_angle)
            except (ValueError, IndexError):
                print(f"Error: Failed to extract numbers from file: {filename}. Check the 'lens sie' line format.")
                return None

    print(f"Error: 'lens sie' line not found after last 'chi^2' in file: {filename}.")
    return None

def rename_file(file_path):
    parts = file_path.split("/")
    lens_system = parts[0]  # Now correctly fetching the first part of the file path
    file_name = parts[-1]

    key = file_name.split("_optresult")[0][3:]
    profile_map = {"S": "SIE", "P": "POW", "N": "NFW"}
    profile = profile_map.get(key[0], "Unknown")

    constraint = "Unknown constraint"
    if key[1:3] == "PF":
        constraint = "QSO Pos + FR"
    elif key[1] == "P":
        constraint = "QSO Pos"

    extras = []
    if "G" in key:
        extras.append("G2")
    if "R" in key:
        extras.append("Shear")

    # Assemble profile descriptions with combinations
    if len(extras) == 2:
        profile_description = f"{profile} + {extras[0]} + {extras[1]}"
    elif len(extras) == 1:
        profile_description = f"{profile} + {extras[0]}"
    else:
        profile_description = profile

    new_name = f"{lens_system} & {profile_description} & {constraint}"
    return new_name

def process_multiple_files(file_list):
    results = []
    for file_name in file_list:
        model_description = rename_file(file_name)
        sie_parameters = extract_sie_parameters(file_name)
        if sie_parameters:
            lens_redshift, velocity_dispersion, x_coord, y_coord, ellipticity, position_angle = sie_parameters
            result_text = f"{model_description} & {lens_redshift:.3f} & {velocity_dispersion:.3f} & {x_coord:.3f} & {y_coord:.3f} & {ellipticity:.3f} & {position_angle:.3f} \\\\"
            results.append(result_text)
            print(f"Model: {file_name}")
            print(f"Lens Redshift: {lens_redshift:.3f}")
            print(f"Velocity Dispersion: {velocity_dispersion:.3f}")
            print(f"X Coordinate: {x_coord:}")
            print(f"Y Coordinate: {y_coord:}")
            print(f"Ellipticity: {ellipticity:.3f}")
            print(f"Position Angle: {position_angle:.3f}\n")
        else:
            results.append(f"{model_description} Failed to extract SIE parameters.")
            print(f"Model: {file_name} Failed to extract SIE parameters.\n")
    return results

# List of file paths to process
file_list = [
    "RXJ0911/SPC/outSP_optresult.dat",
    "RXJ0911/SPC/outSPR_optresult.dat",
    "RXJ0911/SPC/outSPG_optresult.dat",
    "RXJ0911/SPC/outSPGR_optresult.dat",
    "RXJ0911/SPFC/outSPF_optresult.dat",
    "RXJ0911/SPFC/outSPFR_optresult.dat",
    "RXJ0911/SPFC/outSPFG_optresult.dat",
    "RXJ0911/SPFC/outSPFGR_optresult.dat",
    "PSJ1606/SPC/outSP_optresult.dat",
    "PSJ1606/SPC/outSPR_optresult.dat",
    "PSJ1606/SPC/outSPG_optresult.dat",
    "PSJ1606/SPC/outSPGR_optresult.dat",
    "PSJ1606/SPFC/outSPF_optresult.dat",
    "PSJ1606/SPFC/outSPFR_optresult.dat",
    "PSJ1606/SPFC/outSPFG_optresult.dat",
    "PSJ1606/SPFC/outSPFGR_optresult.dat",
    "WFI2033/SP/outSP_optresult.dat",
    "WFI2033/SPR/outSPR_optresult.dat",
    "WFI2033/SPG/outSPG_optresult.dat",
    "WFI2033/SPGR/outSPGR_optresult.dat",
    "WFI2033/SPF/outSPF_optresult.dat",
    "WFI2033/SPFR/outSPFR_optresult.dat",
    "WFI2033/SPFG/outSPFG_optresult.dat",
    "WFI2033/SPFGR/outSPFGR_optresult.dat",
    "SDSSJ1330/SPC/outSP_optresult.dat",
    "SDSSJ1330/SPC/outSPR_optresult.dat",
    "SDSSJ1330/SPFC/outSPF_optresult.dat",
    "SDSSJ1330/SPFC/outSPFR_optresult.dat",
    "WFI2026/SPC/outSP_optresult.dat",
    "WFI2026/SPC/outSPR_optresult.dat",
    "WFI2026/SPFC/outSPF_optresult.dat",
    "WFI2026/SPFC/outSPFR_optresult.dat",
    "WGDJ0405/SPC/outSP_optresult.dat",
    "WGDJ0405/SPC/outSPR_optresult.dat",
    "WGDJ0405/SPFC/outSPF_optresult.dat",
    "WGDJ0405/SPFC/outSPFR_optresult.dat",
    "WGD2038/SPC/outSP_optresult.dat",
    "WGD2038/SPC/outSPR_optresult.dat",
    "WGD2038/SPFC/outSPF_optresult.dat",
    "WGD2038/SPFC/outSPFR_optresult.dat",
]

# Execute the function and print LaTeX-formatted results
latex_results = process_multiple_files(file_list)
for result in latex_results:
    print(result)

