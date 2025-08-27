#!/usr/bin/env python
# coding: utf-8

# In[43]:


# Fixed Combined lens parameter table code
# RXJ0911

def extract_lens_parameters(filename, lens_type):
    with open(filename, 'r') as file:
        lines = file.readlines()

    lens_key_map = {'s': 'sie', 'p': 'pow', 'n': 'anfw'}
    lens_key = lens_key_map[lens_type.lower()]
    shear_present = 'R' in filename.split('/')[-1]
    g2_present = "G" in filename.split('/')[-1]  # Check if 'G' is in the file name

    last_chi2_index = None
    for index, line in enumerate(lines):
        if "chi^2" in line:
            last_chi2_index = index

    if last_chi2_index is None:
        print(f"Error: 'chi^2' line not found in file: {filename}.")
        return None

    result = {}
    lens_found_count = 0
    first_lens_parameters = {}
    g2_parameters = {}
    found_sie_for_g2 = False

    for line in lines[last_chi2_index:]:
        if "lens   " + lens_key in line:
            parts = line.split()
            try:
                current_params = {
                    'x_coord': float(parts[4]),
                    'y_coord': float(parts[5]),
                    'ellipticity': float(parts[6]),
                    'position_angle': float(parts[7]),
                    'velocity_dispersion': float(parts[3]) if lens_type.lower() == 's' else "-",
                    'einstein_radius': float(parts[8]) if lens_type.lower() == 'p' else "-",
                    'pwi': float(parts[9]) if lens_type.lower() == 'p' else "-",
                    'dm_mass': float(parts[3]) if lens_type.lower() == 'n' else "-",
                    'concentration': float(parts[8]) if lens_type.lower() == 'n' else "-"
                }
                lens_found_count += 1
                if lens_found_count == 1:
                    first_lens_parameters = current_params
            except (ValueError, IndexError):
                print(f"Error: Failed to extract parameters from line: {line} in file: {filename}.")
                return None

        # Check specifically for 'lens   sie' after chi^2 for G2 detection in POW and ANFW models
        if "lens   sie" in line and g2_present:
            parts = line.split()
            try:
                g2_parameters = {
                    'x_coord': float(parts[4]),
                    'y_coord': float(parts[5]),
                    'ellipticity': float(parts[6]),
                    'position_angle': float(parts[7]),
                    'velocity_dispersion': float(parts[3])
                }
                found_sie_for_g2 = True
            except (ValueError, IndexError):
                print(f"Error: Failed to extract G2 parameters from line: {line} in file: {filename}.")
                return None

        if shear_present and "lens   pert" in line:
            parts = line.split()
            try:
                result.update({
                    'shear_x_coord': float(parts[4]),
                    'shear_y_coord': float(parts[5]),
                    'shear_position_angle': float(parts[7]),
                    'shear_strength': float(parts[6]),
                    'convergence': float(parts[9])
                })
            except (ValueError, IndexError):
                print(f"Error: Failed to extract shear parameters from line: {line} in file: {filename}.")
                return None

    result.update(first_lens_parameters)  # Add first lens parameters to result

    # Include G2 parameters in the result if applicable and found
    if g2_present and found_sie_for_g2:
        result.update({
            'g2_x_coord': g2_parameters.get('x_coord'),
            'g2_y_coord': g2_parameters.get('y_coord'),
            'g2_ellipticity': g2_parameters.get('ellipticity'),
            'g2_position_angle': g2_parameters.get('position_angle'),
            'g2_velocity_dispersion': g2_parameters.get('velocity_dispersion')
        })

    if not result:
        print(f"Error: '{lens_key}' line not found after last 'chi^2' in file: {filename}.")
        return None
    return result

def get_lens_type_constraint(file_name):
    """
    Extracts the lens type and constraint from the filename.
    """
    base_name = file_name.split('/')[-1]  # Get the file name part
    lens_type = base_name[3]  # 4th character is lens type
    constraint = base_name[4:6]  # 5th and 6th characters are constraint
    return lens_type.lower(), constraint

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
        lens_type, _ = get_lens_type_constraint(file_name)
        model_description = rename_file(file_name)
        parameters = extract_lens_parameters(file_name, lens_type)
        if parameters:
            # Print detailed parameters
            print(f"Model: {file_name}")
            print(f"X Coordinate: {parameters.get('x_coord', '-')}")
            print(f"Y Coordinate: {parameters.get('y_coord', '-')}")
            print(f"Ellipticity: {parameters.get('ellipticity', '-')}")
            print(f"Position Angle: {parameters.get('position_angle', '-')}")
            print(f"Velocity Dispersion: {parameters.get('velocity_dispersion', '-')}")
            print(f"Einstein Radius: {parameters.get('einstein_radius', '-')}")
            print(f"PWI: {parameters.get('pwi', '-')}")
            # Ensure DM Mass is presented in exponential form if it's a float
            dm_mass_formatted = f"{parameters.get('dm_mass', '-')}"
            if isinstance(parameters.get('dm_mass'), float):
                dm_mass_formatted = f"{parameters['dm_mass']:e}"
            print(f"DM Mass: {dm_mass_formatted}")
            print(f"Concentration Parameter: {parameters.get('concentration', '-')}")
            if 'shear_x_coord' in parameters:
                print(f"Shear X Coordinate: {parameters['shear_x_coord']}")
                print(f"Shear Y Coordinate: {parameters['shear_y_coord']}")
                print(f"Shear Position Angle: {parameters['shear_position_angle']}")
                print(f"Shear Strength: {parameters['shear_strength']}")
                print(f"Convergence: {parameters['convergence']}")
            if 'g2_x_coord' in parameters:
                print(f"G2 X Coordinate: {parameters['g2_x_coord']}")
                print(f"G2 Y Coordinate: {parameters['g2_y_coord']}")
                print(f"G2 Ellipticity: {parameters['g2_ellipticity']}")
                print(f"G2 Position Angle: {parameters['g2_position_angle']}")
                print(f"G2 Velocity Dispersion: {parameters['g2_velocity_dispersion']}")
            print()

            # Continue with LaTeX output formatting
            result_text = f"{model_description} & {parameters.get('x_coord', '-'):0.3f} & {parameters.get('y_coord', '-'):0.3f} & {parameters.get('ellipticity', '-'):0.3f} & {parameters.get('position_angle', '-'):0.3f} & {parameters.get('velocity_dispersion', '-')} & {parameters.get('einstein_radius', '-')} & {parameters.get('pwi', '-')} & {dm_mass_formatted} & {parameters.get('concentration', '-')} & - & - \\\\"
            if 'shear_x_coord' in parameters:
                result_text += f"\n - & Shear & - & {parameters['shear_x_coord']:0.3f} & {parameters['shear_y_coord']:0.3f} & - & {parameters['shear_position_angle']:0.3f} & - & - & - & - & - & {parameters['shear_strength']:0.3f} & {parameters['convergence']:0.3f} \\\\"
            if 'g2_x_coord' in parameters:
                result_text += f"\n - & G2 & - & {parameters['g2_x_coord']:0.3f} & {parameters['g2_y_coord']:0.3f} & {parameters['g2_ellipticity']:0.3f} & {parameters['g2_position_angle']:0.3f} & {parameters['g2_velocity_dispersion']} & - & - & - & - & - & - \\\\"
            result_text += "\n\\hline"
            results.append(result_text)
        else:
            result_text = f"{model_description} & No valid data for {lens_type.upper()} \\\\"
            result_text += "\n\\hline"
            results.append(result_text)
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
    "RXJ0911/NPFC/outNPFGR_optresult.dat",
]

# Execute the function and print LaTeX-formatted results
results = process_multiple_files(file_list)
for result in results:
    print(result)


# In[44]:


# Fixed Combined lens parameter table code
# PSJ1606

def extract_lens_parameters(filename, lens_type):
    with open(filename, 'r') as file:
        lines = file.readlines()

    lens_key_map = {'s': 'sie', 'p': 'pow', 'n': 'anfw'}
    lens_key = lens_key_map[lens_type.lower()]
    shear_present = 'R' in filename.split('/')[-1]
    g2_present = "G" in filename.split('/')[-1]  # Check if 'G' is in the file name

    last_chi2_index = None
    for index, line in enumerate(lines):
        if "chi^2" in line:
            last_chi2_index = index

    if last_chi2_index is None:
        print(f"Error: 'chi^2' line not found in file: {filename}.")
        return None

    result = {}
    lens_found_count = 0
    first_lens_parameters = {}
    g2_parameters = {}
    found_sie_for_g2 = False

    for line in lines[last_chi2_index:]:
        if "lens   " + lens_key in line:
            parts = line.split()
            try:
                current_params = {
                    'x_coord': float(parts[4]),
                    'y_coord': float(parts[5]),
                    'ellipticity': float(parts[6]),
                    'position_angle': float(parts[7]),
                    'velocity_dispersion': float(parts[3]) if lens_type.lower() == 's' else "-",
                    'einstein_radius': float(parts[8]) if lens_type.lower() == 'p' else "-",
                    'pwi': float(parts[9]) if lens_type.lower() == 'p' else "-",
                    'dm_mass': float(parts[3]) if lens_type.lower() == 'n' else "-",
                    'concentration': float(parts[8]) if lens_type.lower() == 'n' else "-"
                }
                lens_found_count += 1
                if lens_found_count == 1:
                    first_lens_parameters = current_params
            except (ValueError, IndexError):
                print(f"Error: Failed to extract parameters from line: {line} in file: {filename}.")
                return None

        # Check specifically for 'lens   sie' after chi^2 for G2 detection in POW and ANFW models
        if "lens   sie" in line and g2_present:
            parts = line.split()
            try:
                g2_parameters = {
                    'x_coord': float(parts[4]),
                    'y_coord': float(parts[5]),
                    'ellipticity': float(parts[6]),
                    'position_angle': float(parts[7]),
                    'velocity_dispersion': float(parts[3])
                }
                found_sie_for_g2 = True
            except (ValueError, IndexError):
                print(f"Error: Failed to extract G2 parameters from line: {line} in file: {filename}.")
                return None

        if shear_present and "lens   pert" in line:
            parts = line.split()
            try:
                result.update({
                    'shear_x_coord': float(parts[4]),
                    'shear_y_coord': float(parts[5]),
                    'shear_position_angle': float(parts[7]),
                    'shear_strength': float(parts[6]),
                    'convergence': float(parts[9])
                })
            except (ValueError, IndexError):
                print(f"Error: Failed to extract shear parameters from line: {line} in file: {filename}.")
                return None

    result.update(first_lens_parameters)  # Add first lens parameters to result

    # Include G2 parameters in the result if applicable and found
    if g2_present and found_sie_for_g2:
        result.update({
            'g2_x_coord': g2_parameters.get('x_coord'),
            'g2_y_coord': g2_parameters.get('y_coord'),
            'g2_ellipticity': g2_parameters.get('ellipticity'),
            'g2_position_angle': g2_parameters.get('position_angle'),
            'g2_velocity_dispersion': g2_parameters.get('velocity_dispersion')
        })

    if not result:
        print(f"Error: '{lens_key}' line not found after last 'chi^2' in file: {filename}.")
        return None
    return result

def get_lens_type_constraint(file_name):
    """
    Extracts the lens type and constraint from the filename.
    """
    base_name = file_name.split('/')[-1]  # Get the file name part
    lens_type = base_name[3]  # 4th character is lens type
    constraint = base_name[4:6]  # 5th and 6th characters are constraint
    return lens_type.lower(), constraint

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
        lens_type, _ = get_lens_type_constraint(file_name)
        model_description = rename_file(file_name)
        parameters = extract_lens_parameters(file_name, lens_type)
        if parameters:
            # Print detailed parameters
            print(f"Model: {file_name}")
            print(f"X Coordinate: {parameters.get('x_coord', '-')}")
            print(f"Y Coordinate: {parameters.get('y_coord', '-')}")
            print(f"Ellipticity: {parameters.get('ellipticity', '-')}")
            print(f"Position Angle: {parameters.get('position_angle', '-')}")
            print(f"Velocity Dispersion: {parameters.get('velocity_dispersion', '-')}")
            print(f"Einstein Radius: {parameters.get('einstein_radius', '-')}")
            print(f"PWI: {parameters.get('pwi', '-')}")
            # Ensure DM Mass is presented in exponential form if it's a float
            dm_mass_formatted = f"{parameters.get('dm_mass', '-')}"
            if isinstance(parameters.get('dm_mass'), float):
                dm_mass_formatted = f"{parameters['dm_mass']:e}"
            print(f"DM Mass: {dm_mass_formatted}")
            print(f"Concentration Parameter: {parameters.get('concentration', '-')}")
            if 'shear_x_coord' in parameters:
                print(f"Shear X Coordinate: {parameters['shear_x_coord']}")
                print(f"Shear Y Coordinate: {parameters['shear_y_coord']}")
                print(f"Shear Position Angle: {parameters['shear_position_angle']}")
                print(f"Shear Strength: {parameters['shear_strength']}")
                print(f"Convergence: {parameters['convergence']}")
            if 'g2_x_coord' in parameters:
                print(f"G2 X Coordinate: {parameters['g2_x_coord']}")
                print(f"G2 Y Coordinate: {parameters['g2_y_coord']}")
                print(f"G2 Ellipticity: {parameters['g2_ellipticity']}")
                print(f"G2 Position Angle: {parameters['g2_position_angle']}")
                print(f"G2 Velocity Dispersion: {parameters['g2_velocity_dispersion']}")
            print()

            # Continue with LaTeX output formatting
            result_text = f"{model_description} & {parameters.get('x_coord', '-'):0.3f} & {parameters.get('y_coord', '-'):0.3f} & {parameters.get('ellipticity', '-'):0.3f} & {parameters.get('position_angle', '-'):0.3f} & {parameters.get('velocity_dispersion', '-')} & {parameters.get('einstein_radius', '-')} & {parameters.get('pwi', '-')} & {dm_mass_formatted} & {parameters.get('concentration', '-')} & - & - \\\\"
            if 'shear_x_coord' in parameters:
                result_text += f"\n - & Shear & - & {parameters['shear_x_coord']:0.3f} & {parameters['shear_y_coord']:0.3f} & - & {parameters['shear_position_angle']:0.3f} & - & - & - & - & - & {parameters['shear_strength']:0.3f} & {parameters['convergence']:0.3f} \\\\"
            if 'g2_x_coord' in parameters:
                result_text += f"\n - & G2 & - & {parameters['g2_x_coord']:0.3f} & {parameters['g2_y_coord']:0.3f} & {parameters['g2_ellipticity']:0.3f} & {parameters['g2_position_angle']:0.3f} & {parameters['g2_velocity_dispersion']} & - & - & - & - & - & - \\\\"
            result_text += "\n\\hline"
            results.append(result_text)
        else:
            result_text = f"{model_description} & No valid data for {lens_type.upper()} \\\\"
            result_text += "\n\\hline"
            results.append(result_text)
    return results

# List of file paths to process
file_list = [
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

# Execute the function and print LaTeX-formatted results
results = process_multiple_files(file_list)
for result in results:
    print(result)


# In[45]:


# Fixed Combined lens parameter table code
# WFI2033

def extract_lens_parameters(filename, lens_type):
    with open(filename, 'r') as file:
        lines = file.readlines()

    lens_key_map = {'s': 'sie', 'p': 'pow', 'n': 'anfw'}
    lens_key = lens_key_map[lens_type.lower()]
    shear_present = 'R' in filename.split('/')[-1]
    g2_present = "G" in filename.split('/')[-1]  # Check if 'G' is in the file name

    last_chi2_index = None
    for index, line in enumerate(lines):
        if "chi^2" in line:
            last_chi2_index = index

    if last_chi2_index is None:
        print(f"Error: 'chi^2' line not found in file: {filename}.")
        return None

    result = {}
    lens_found_count = 0
    first_lens_parameters = {}
    g2_parameters = {}
    found_sie_for_g2 = False

    for line in lines[last_chi2_index:]:
        if "lens   " + lens_key in line:
            parts = line.split()
            try:
                current_params = {
                    'x_coord': float(parts[4]),
                    'y_coord': float(parts[5]),
                    'ellipticity': float(parts[6]),
                    'position_angle': float(parts[7]),
                    'velocity_dispersion': float(parts[3]) if lens_type.lower() == 's' else "-",
                    'einstein_radius': float(parts[8]) if lens_type.lower() == 'p' else "-",
                    'pwi': float(parts[9]) if lens_type.lower() == 'p' else "-",
                    'dm_mass': float(parts[3]) if lens_type.lower() == 'n' else "-",
                    'concentration': float(parts[8]) if lens_type.lower() == 'n' else "-"
                }
                lens_found_count += 1
                if lens_found_count == 1:
                    first_lens_parameters = current_params
            except (ValueError, IndexError):
                print(f"Error: Failed to extract parameters from line: {line} in file: {filename}.")
                return None

        # Check specifically for 'lens   sie' after chi^2 for G2 detection in POW and ANFW models
        if "lens   sie" in line and g2_present:
            parts = line.split()
            try:
                g2_parameters = {
                    'x_coord': float(parts[4]),
                    'y_coord': float(parts[5]),
                    'ellipticity': float(parts[6]),
                    'position_angle': float(parts[7]),
                    'velocity_dispersion': float(parts[3])
                }
                found_sie_for_g2 = True
            except (ValueError, IndexError):
                print(f"Error: Failed to extract G2 parameters from line: {line} in file: {filename}.")
                return None

        if shear_present and "lens   pert" in line:
            parts = line.split()
            try:
                result.update({
                    'shear_x_coord': float(parts[4]),
                    'shear_y_coord': float(parts[5]),
                    'shear_position_angle': float(parts[7]),
                    'shear_strength': float(parts[6]),
                    'convergence': float(parts[9])
                })
            except (ValueError, IndexError):
                print(f"Error: Failed to extract shear parameters from line: {line} in file: {filename}.")
                return None

    result.update(first_lens_parameters)  # Add first lens parameters to result

    # Include G2 parameters in the result if applicable and found
    if g2_present and found_sie_for_g2:
        result.update({
            'g2_x_coord': g2_parameters.get('x_coord'),
            'g2_y_coord': g2_parameters.get('y_coord'),
            'g2_ellipticity': g2_parameters.get('ellipticity'),
            'g2_position_angle': g2_parameters.get('position_angle'),
            'g2_velocity_dispersion': g2_parameters.get('velocity_dispersion')
        })

    if not result:
        print(f"Error: '{lens_key}' line not found after last 'chi^2' in file: {filename}.")
        return None
    return result

def get_lens_type_constraint(file_name):
    """
    Extracts the lens type and constraint from the filename.
    """
    base_name = file_name.split('/')[-1]  # Get the file name part
    lens_type = base_name[3]  # 4th character is lens type
    constraint = base_name[4:6]  # 5th and 6th characters are constraint
    return lens_type.lower(), constraint

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
        lens_type, _ = get_lens_type_constraint(file_name)
        model_description = rename_file(file_name)
        parameters = extract_lens_parameters(file_name, lens_type)
        if parameters:
            # Print detailed parameters
            print(f"Model: {file_name}")
            print(f"X Coordinate: {parameters.get('x_coord', '-')}")
            print(f"Y Coordinate: {parameters.get('y_coord', '-')}")
            print(f"Ellipticity: {parameters.get('ellipticity', '-')}")
            print(f"Position Angle: {parameters.get('position_angle', '-')}")
            print(f"Velocity Dispersion: {parameters.get('velocity_dispersion', '-')}")
            print(f"Einstein Radius: {parameters.get('einstein_radius', '-')}")
            print(f"PWI: {parameters.get('pwi', '-')}")
            # Ensure DM Mass is presented in exponential form if it's a float
            dm_mass_formatted = f"{parameters.get('dm_mass', '-')}"
            if isinstance(parameters.get('dm_mass'), float):
                dm_mass_formatted = f"{parameters['dm_mass']:e}"
            print(f"DM Mass: {dm_mass_formatted}")
            print(f"Concentration Parameter: {parameters.get('concentration', '-')}")
            if 'shear_x_coord' in parameters:
                print(f"Shear X Coordinate: {parameters['shear_x_coord']}")
                print(f"Shear Y Coordinate: {parameters['shear_y_coord']}")
                print(f"Shear Position Angle: {parameters['shear_position_angle']}")
                print(f"Shear Strength: {parameters['shear_strength']}")
                print(f"Convergence: {parameters['convergence']}")
            if 'g2_x_coord' in parameters:
                print(f"G2 X Coordinate: {parameters['g2_x_coord']}")
                print(f"G2 Y Coordinate: {parameters['g2_y_coord']}")
                print(f"G2 Ellipticity: {parameters['g2_ellipticity']}")
                print(f"G2 Position Angle: {parameters['g2_position_angle']}")
                print(f"G2 Velocity Dispersion: {parameters['g2_velocity_dispersion']}")
            print()

            # Continue with LaTeX output formatting
            result_text = f"{model_description} & {parameters.get('x_coord', '-'):0.3f} & {parameters.get('y_coord', '-'):0.3f} & {parameters.get('ellipticity', '-'):0.3f} & {parameters.get('position_angle', '-'):0.3f} & {parameters.get('velocity_dispersion', '-')} & {parameters.get('einstein_radius', '-')} & {parameters.get('pwi', '-')} & {dm_mass_formatted} & {parameters.get('concentration', '-')} & - & - \\\\"
            if 'shear_x_coord' in parameters:
                result_text += f"\n - & Shear & - & {parameters['shear_x_coord']:0.3f} & {parameters['shear_y_coord']:0.3f} & - & {parameters['shear_position_angle']:0.3f} & - & - & - & - & - & {parameters['shear_strength']:0.3f} & {parameters['convergence']:0.3f} \\\\"
            if 'g2_x_coord' in parameters:
                result_text += f"\n - & G2 & - & {parameters['g2_x_coord']:0.3f} & {parameters['g2_y_coord']:0.3f} & {parameters['g2_ellipticity']:0.3f} & {parameters['g2_position_angle']:0.3f} & {parameters['g2_velocity_dispersion']} & - & - & - & - & - & - \\\\"
            result_text += "\n\\hline"
            results.append(result_text)
        else:
            result_text = f"{model_description} & No valid data for {lens_type.upper()} \\\\"
            result_text += "\n\\hline"
            results.append(result_text)
    return results

# List of file paths to process
file_list = [
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

# Execute the function and print LaTeX-formatted results
results = process_multiple_files(file_list)
for result in results:
    print(result)


# In[1]:


# Fixed Combined lens parameter table code
# SDSSJ1330

def extract_lens_parameters(filename, lens_type):
    with open(filename, 'r') as file:
        lines = file.readlines()

    lens_key_map = {'s': 'sie', 'p': 'pow', 'n': 'anfw'}
    lens_key = lens_key_map[lens_type.lower()]
    shear_present = 'R' in filename.split('/')[-1]
    g2_present = "G" in filename.split('/')[-1]  # Check if 'G' is in the file name

    last_chi2_index = None
    for index, line in enumerate(lines):
        if "chi^2" in line:
            last_chi2_index = index

    if last_chi2_index is None:
        print(f"Error: 'chi^2' line not found in file: {filename}.")
        return None

    result = {}
    lens_found_count = 0
    first_lens_parameters = {}
    g2_parameters = {}
    found_sie_for_g2 = False

    for line in lines[last_chi2_index:]:
        if "lens   " + lens_key in line:
            parts = line.split()
            try:
                current_params = {
                    'x_coord': float(parts[4]),
                    'y_coord': float(parts[5]),
                    'ellipticity': float(parts[6]),
                    'position_angle': float(parts[7]),
                    'velocity_dispersion': float(parts[3]) if lens_type.lower() == 's' else "-",
                    'einstein_radius': float(parts[8]) if lens_type.lower() == 'p' else "-",
                    'pwi': float(parts[9]) if lens_type.lower() == 'p' else "-",
                    'dm_mass': float(parts[3]) if lens_type.lower() == 'n' else "-",
                    'concentration': float(parts[8]) if lens_type.lower() == 'n' else "-"
                }
                lens_found_count += 1
                if lens_found_count == 1:
                    first_lens_parameters = current_params
            except (ValueError, IndexError):
                print(f"Error: Failed to extract parameters from line: {line} in file: {filename}.")
                return None

        # Check specifically for 'lens   sie' after chi^2 for G2 detection in POW and ANFW models
        if "lens   sie" in line and g2_present:
            parts = line.split()
            try:
                g2_parameters = {
                    'x_coord': float(parts[4]),
                    'y_coord': float(parts[5]),
                    'ellipticity': float(parts[6]),
                    'position_angle': float(parts[7]),
                    'velocity_dispersion': float(parts[3])
                }
                found_sie_for_g2 = True
            except (ValueError, IndexError):
                print(f"Error: Failed to extract G2 parameters from line: {line} in file: {filename}.")
                return None

        if shear_present and "lens   pert" in line:
            parts = line.split()
            try:
                result.update({
                    'shear_x_coord': float(parts[4]),
                    'shear_y_coord': float(parts[5]),
                    'shear_position_angle': float(parts[7]),
                    'shear_strength': float(parts[6]),
                    'convergence': float(parts[9])
                })
            except (ValueError, IndexError):
                print(f"Error: Failed to extract shear parameters from line: {line} in file: {filename}.")
                return None

    result.update(first_lens_parameters)  # Add first lens parameters to result

    # Include G2 parameters in the result if applicable and found
    if g2_present and found_sie_for_g2:
        result.update({
            'g2_x_coord': g2_parameters.get('x_coord'),
            'g2_y_coord': g2_parameters.get('y_coord'),
            'g2_ellipticity': g2_parameters.get('ellipticity'),
            'g2_position_angle': g2_parameters.get('position_angle'),
            'g2_velocity_dispersion': g2_parameters.get('velocity_dispersion')
        })

    if not result:
        print(f"Error: '{lens_key}' line not found after last 'chi^2' in file: {filename}.")
        return None
    return result

def get_lens_type_constraint(file_name):
    """
    Extracts the lens type and constraint from the filename.
    """
    base_name = file_name.split('/')[-1]  # Get the file name part
    lens_type = base_name[3]  # 4th character is lens type
    constraint = base_name[4:6]  # 5th and 6th characters are constraint
    return lens_type.lower(), constraint

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
        lens_type, _ = get_lens_type_constraint(file_name)
        model_description = rename_file(file_name)
        parameters = extract_lens_parameters(file_name, lens_type)
        if parameters:
            # Print detailed parameters
            print(f"Model: {file_name}")
            print(f"X Coordinate: {parameters.get('x_coord', '-')}")
            print(f"Y Coordinate: {parameters.get('y_coord', '-')}")
            print(f"Ellipticity: {parameters.get('ellipticity', '-')}")
            print(f"Position Angle: {parameters.get('position_angle', '-')}")
            print(f"Velocity Dispersion: {parameters.get('velocity_dispersion', '-')}")
            print(f"Einstein Radius: {parameters.get('einstein_radius', '-')}")
            print(f"PWI: {parameters.get('pwi', '-')}")
            # Ensure DM Mass is presented in exponential form if it's a float
            dm_mass_formatted = f"{parameters.get('dm_mass', '-')}"
            if isinstance(parameters.get('dm_mass'), float):
                dm_mass_formatted = f"{parameters['dm_mass']:e}"
            print(f"DM Mass: {dm_mass_formatted}")
            print(f"Concentration Parameter: {parameters.get('concentration', '-')}")
            if 'shear_x_coord' in parameters:
                print(f"Shear X Coordinate: {parameters['shear_x_coord']}")
                print(f"Shear Y Coordinate: {parameters['shear_y_coord']}")
                print(f"Shear Position Angle: {parameters['shear_position_angle']}")
                print(f"Shear Strength: {parameters['shear_strength']}")
                print(f"Convergence: {parameters['convergence']}")
            if 'g2_x_coord' in parameters:
                print(f"G2 X Coordinate: {parameters['g2_x_coord']}")
                print(f"G2 Y Coordinate: {parameters['g2_y_coord']}")
                print(f"G2 Ellipticity: {parameters['g2_ellipticity']}")
                print(f"G2 Position Angle: {parameters['g2_position_angle']}")
                print(f"G2 Velocity Dispersion: {parameters['g2_velocity_dispersion']}")
            print()

            # Continue with LaTeX output formatting
            result_text = f"{model_description} & {parameters.get('x_coord', '-'):0.3f} & {parameters.get('y_coord', '-'):0.3f} & {parameters.get('ellipticity', '-'):0.3f} & {parameters.get('position_angle', '-'):0.3f} & {parameters.get('velocity_dispersion', '-')} & {parameters.get('einstein_radius', '-')} & {parameters.get('pwi', '-')} & {dm_mass_formatted} & {parameters.get('concentration', '-')} & - & - \\\\"
            if 'shear_x_coord' in parameters:
                result_text += f"\n - & Shear & - & {parameters['shear_x_coord']:0.3f} & {parameters['shear_y_coord']:0.3f} & - & {parameters['shear_position_angle']:0.3f} & - & - & - & - & - & {parameters['shear_strength']:0.3f} & {parameters['convergence']:0.3f} \\\\"
            if 'g2_x_coord' in parameters:
                result_text += f"\n - & G2 & - & {parameters['g2_x_coord']:0.3f} & {parameters['g2_y_coord']:0.3f} & {parameters['g2_ellipticity']:0.3f} & {parameters['g2_position_angle']:0.3f} & {parameters['g2_velocity_dispersion']} & - & - & - & - & - & - \\\\"
            result_text += "\n\\hline"
            results.append(result_text)
        else:
            result_text = f"{model_description} & No valid data for {lens_type.upper()} \\\\"
            result_text += "\n\\hline"
            results.append(result_text)
    return results

# List of file paths to process
file_list = [
#    "SDSSJ1330/SPC/outSP_optresult.dat",
#    "SDSSJ1330/SPC/outSPR_optresult.dat",
#    "SDSSJ1330/SPFC/outSPF_optresult.dat",
#    "SDSSJ1330/SPFC/outSPFR_optresult.dat",
#    "SDSSJ1330/PPC/outPP_optresult.dat",
    "SDSSJ1330/PPC/outPP_N_optresult.dat",
#    "SDSSJ1330/PPC/outPPR_optresult.dat",
    "SDSSJ1330/PPC/outPPR_N_optresult.dat",
#    "SDSSJ1330/PPFC/outPPF_optresult.dat",
#    "SDSSJ1330/PPFC/outPPFR_optresult.dat",
#    "SDSSJ1330/NPC/outNP_optresult.dat",
#    "SDSSJ1330/NPC/outNPR_optresult.dat",
#    "SDSSJ1330/NPFC/outNPF_optresult.dat",
#    "SDSSJ1330/NPFC/outNPFR_optresult.dat",
]

# Execute the function and print LaTeX-formatted results
results = process_multiple_files(file_list)
for result in results:
    print(result)


# In[2]:


# Fixed Combined lens parameter table code
# WFI2026

def extract_lens_parameters(filename, lens_type):
    with open(filename, 'r') as file:
        lines = file.readlines()

    lens_key_map = {'s': 'sie', 'p': 'pow', 'n': 'anfw'}
    lens_key = lens_key_map[lens_type.lower()]
    shear_present = 'R' in filename.split('/')[-1]
    g2_present = "G" in filename.split('/')[-1]  # Check if 'G' is in the file name

    last_chi2_index = None
    for index, line in enumerate(lines):
        if "chi^2" in line:
            last_chi2_index = index

    if last_chi2_index is None:
        print(f"Error: 'chi^2' line not found in file: {filename}.")
        return None

    result = {}
    lens_found_count = 0
    first_lens_parameters = {}
    g2_parameters = {}
    found_sie_for_g2 = False

    for line in lines[last_chi2_index:]:
        if "lens   " + lens_key in line:
            parts = line.split()
            try:
                current_params = {
                    'x_coord': float(parts[4]),
                    'y_coord': float(parts[5]),
                    'ellipticity': float(parts[6]),
                    'position_angle': float(parts[7]),
                    'velocity_dispersion': float(parts[3]) if lens_type.lower() == 's' else "-",
                    'einstein_radius': float(parts[8]) if lens_type.lower() == 'p' else "-",
                    'pwi': float(parts[9]) if lens_type.lower() == 'p' else "-",
                    'dm_mass': float(parts[3]) if lens_type.lower() == 'n' else "-",
                    'concentration': float(parts[8]) if lens_type.lower() == 'n' else "-"
                }
                lens_found_count += 1
                if lens_found_count == 1:
                    first_lens_parameters = current_params
            except (ValueError, IndexError):
                print(f"Error: Failed to extract parameters from line: {line} in file: {filename}.")
                return None

        # Check specifically for 'lens   sie' after chi^2 for G2 detection in POW and ANFW models
        if "lens   sie" in line and g2_present:
            parts = line.split()
            try:
                g2_parameters = {
                    'x_coord': float(parts[4]),
                    'y_coord': float(parts[5]),
                    'ellipticity': float(parts[6]),
                    'position_angle': float(parts[7]),
                    'velocity_dispersion': float(parts[3])
                }
                found_sie_for_g2 = True
            except (ValueError, IndexError):
                print(f"Error: Failed to extract G2 parameters from line: {line} in file: {filename}.")
                return None

        if shear_present and "lens   pert" in line:
            parts = line.split()
            try:
                result.update({
                    'shear_x_coord': float(parts[4]),
                    'shear_y_coord': float(parts[5]),
                    'shear_position_angle': float(parts[7]),
                    'shear_strength': float(parts[6]),
                    'convergence': float(parts[9])
                })
            except (ValueError, IndexError):
                print(f"Error: Failed to extract shear parameters from line: {line} in file: {filename}.")
                return None

    result.update(first_lens_parameters)  # Add first lens parameters to result

    # Include G2 parameters in the result if applicable and found
    if g2_present and found_sie_for_g2:
        result.update({
            'g2_x_coord': g2_parameters.get('x_coord'),
            'g2_y_coord': g2_parameters.get('y_coord'),
            'g2_ellipticity': g2_parameters.get('ellipticity'),
            'g2_position_angle': g2_parameters.get('position_angle'),
            'g2_velocity_dispersion': g2_parameters.get('velocity_dispersion')
        })

    if not result:
        print(f"Error: '{lens_key}' line not found after last 'chi^2' in file: {filename}.")
        return None
    return result

def get_lens_type_constraint(file_name):
    """
    Extracts the lens type and constraint from the filename.
    """
    base_name = file_name.split('/')[-1]  # Get the file name part
    lens_type = base_name[3]  # 4th character is lens type
    constraint = base_name[4:6]  # 5th and 6th characters are constraint
    return lens_type.lower(), constraint

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
        lens_type, _ = get_lens_type_constraint(file_name)
        model_description = rename_file(file_name)
        parameters = extract_lens_parameters(file_name, lens_type)
        if parameters:
            # Print detailed parameters
            print(f"Model: {file_name}")
            print(f"X Coordinate: {parameters.get('x_coord', '-')}")
            print(f"Y Coordinate: {parameters.get('y_coord', '-')}")
            print(f"Ellipticity: {parameters.get('ellipticity', '-')}")
            print(f"Position Angle: {parameters.get('position_angle', '-')}")
            print(f"Velocity Dispersion: {parameters.get('velocity_dispersion', '-')}")
            print(f"Einstein Radius: {parameters.get('einstein_radius', '-')}")
            print(f"PWI: {parameters.get('pwi', '-')}")
            # Ensure DM Mass is presented in exponential form if it's a float
            dm_mass_formatted = f"{parameters.get('dm_mass', '-')}"
            if isinstance(parameters.get('dm_mass'), float):
                dm_mass_formatted = f"{parameters['dm_mass']:e}"
            print(f"DM Mass: {dm_mass_formatted}")
            print(f"Concentration Parameter: {parameters.get('concentration', '-')}")
            if 'shear_x_coord' in parameters:
                print(f"Shear X Coordinate: {parameters['shear_x_coord']}")
                print(f"Shear Y Coordinate: {parameters['shear_y_coord']}")
                print(f"Shear Position Angle: {parameters['shear_position_angle']}")
                print(f"Shear Strength: {parameters['shear_strength']}")
                print(f"Convergence: {parameters['convergence']}")
            if 'g2_x_coord' in parameters:
                print(f"G2 X Coordinate: {parameters['g2_x_coord']}")
                print(f"G2 Y Coordinate: {parameters['g2_y_coord']}")
                print(f"G2 Ellipticity: {parameters['g2_ellipticity']}")
                print(f"G2 Position Angle: {parameters['g2_position_angle']}")
                print(f"G2 Velocity Dispersion: {parameters['g2_velocity_dispersion']}")
            print()

            # Continue with LaTeX output formatting
            result_text = f"{model_description} & {parameters.get('x_coord', '-'):0.3f} & {parameters.get('y_coord', '-'):0.3f} & {parameters.get('ellipticity', '-'):0.3f} & {parameters.get('position_angle', '-'):0.3f} & {parameters.get('velocity_dispersion', '-')} & {parameters.get('einstein_radius', '-')} & {parameters.get('pwi', '-')} & {dm_mass_formatted} & {parameters.get('concentration', '-')} & - & - \\\\"
            if 'shear_x_coord' in parameters:
                result_text += f"\n - & Shear & - & {parameters['shear_x_coord']:0.3f} & {parameters['shear_y_coord']:0.3f} & - & {parameters['shear_position_angle']:0.3f} & - & - & - & - & - & {parameters['shear_strength']:0.3f} & {parameters['convergence']:0.3f} \\\\"
            if 'g2_x_coord' in parameters:
                result_text += f"\n - & G2 & - & {parameters['g2_x_coord']:0.3f} & {parameters['g2_y_coord']:0.3f} & {parameters['g2_ellipticity']:0.3f} & {parameters['g2_position_angle']:0.3f} & {parameters['g2_velocity_dispersion']} & - & - & - & - & - & - \\\\"
            result_text += "\n\\hline"
            results.append(result_text)
        else:
            result_text = f"{model_description} & No valid data for {lens_type.upper()} \\\\"
            result_text += "\n\\hline"
            results.append(result_text)
    return results

# List of file paths to process
file_list = [
#    "WFI2026/SPC/outSP_optresult.dat",
#    "WFI2026/SPC/outSPR_optresult.dat",
#    "WFI2026/SPFC/outSPF_optresult.dat",
#    "WFI2026/SPFC/outSPFR_optresult.dat",
#    "WFI2026/PPC/outPP_optresult.dat",
    "WFI2026/PPC/outPP_N_optresult.dat",
#    "WFI2026/PPC/outPPR_optresult.dat",
    "WFI2026/PPC/outPPR_N_optresult.dat",
#    "WFI2026/PPFC/outPPF_optresult.dat",
#    "WFI2026/PPFC/outPPFR_optresult.dat",
#    "WFI2026/NPC/outNP_optresult.dat",
#    "WFI2026/NPC/outNPR_optresult.dat",
#    "WFI2026/NPFC/outNPF_optresult.dat",
#    "WFI2026/NPFC/outNPFR_optresult.dat",
]

# Execute the function and print LaTeX-formatted results
results = process_multiple_files(file_list)
for result in results:
    print(result)


# In[3]:


# Fixed Combined lens parameter table code
# WGDJ0405

def extract_lens_parameters(filename, lens_type):
    with open(filename, 'r') as file:
        lines = file.readlines()

    lens_key_map = {'s': 'sie', 'p': 'pow', 'n': 'anfw'}
    lens_key = lens_key_map[lens_type.lower()]
    shear_present = 'R' in filename.split('/')[-1]
    g2_present = "G" in filename.split('/')[-1]  # Check if 'G' is in the file name

    last_chi2_index = None
    for index, line in enumerate(lines):
        if "chi^2" in line:
            last_chi2_index = index

    if last_chi2_index is None:
        print(f"Error: 'chi^2' line not found in file: {filename}.")
        return None

    result = {}
    lens_found_count = 0
    first_lens_parameters = {}
    g2_parameters = {}
    found_sie_for_g2 = False

    for line in lines[last_chi2_index:]:
        if "lens   " + lens_key in line:
            parts = line.split()
            try:
                current_params = {
                    'x_coord': float(parts[4]),
                    'y_coord': float(parts[5]),
                    'ellipticity': float(parts[6]),
                    'position_angle': float(parts[7]),
                    'velocity_dispersion': float(parts[3]) if lens_type.lower() == 's' else "-",
                    'einstein_radius': float(parts[8]) if lens_type.lower() == 'p' else "-",
                    'pwi': float(parts[9]) if lens_type.lower() == 'p' else "-",
                    'dm_mass': float(parts[3]) if lens_type.lower() == 'n' else "-",
                    'concentration': float(parts[8]) if lens_type.lower() == 'n' else "-"
                }
                lens_found_count += 1
                if lens_found_count == 1:
                    first_lens_parameters = current_params
            except (ValueError, IndexError):
                print(f"Error: Failed to extract parameters from line: {line} in file: {filename}.")
                return None

        # Check specifically for 'lens   sie' after chi^2 for G2 detection in POW and ANFW models
        if "lens   sie" in line and g2_present:
            parts = line.split()
            try:
                g2_parameters = {
                    'x_coord': float(parts[4]),
                    'y_coord': float(parts[5]),
                    'ellipticity': float(parts[6]),
                    'position_angle': float(parts[7]),
                    'velocity_dispersion': float(parts[3])
                }
                found_sie_for_g2 = True
            except (ValueError, IndexError):
                print(f"Error: Failed to extract G2 parameters from line: {line} in file: {filename}.")
                return None

        if shear_present and "lens   pert" in line:
            parts = line.split()
            try:
                result.update({
                    'shear_x_coord': float(parts[4]),
                    'shear_y_coord': float(parts[5]),
                    'shear_position_angle': float(parts[7]),
                    'shear_strength': float(parts[6]),
                    'convergence': float(parts[9])
                })
            except (ValueError, IndexError):
                print(f"Error: Failed to extract shear parameters from line: {line} in file: {filename}.")
                return None

    result.update(first_lens_parameters)  # Add first lens parameters to result

    # Include G2 parameters in the result if applicable and found
    if g2_present and found_sie_for_g2:
        result.update({
            'g2_x_coord': g2_parameters.get('x_coord'),
            'g2_y_coord': g2_parameters.get('y_coord'),
            'g2_ellipticity': g2_parameters.get('ellipticity'),
            'g2_position_angle': g2_parameters.get('position_angle'),
            'g2_velocity_dispersion': g2_parameters.get('velocity_dispersion')
        })

    if not result:
        print(f"Error: '{lens_key}' line not found after last 'chi^2' in file: {filename}.")
        return None
    return result

def get_lens_type_constraint(file_name):
    """
    Extracts the lens type and constraint from the filename.
    """
    base_name = file_name.split('/')[-1]  # Get the file name part
    lens_type = base_name[3]  # 4th character is lens type
    constraint = base_name[4:6]  # 5th and 6th characters are constraint
    return lens_type.lower(), constraint

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
        lens_type, _ = get_lens_type_constraint(file_name)
        model_description = rename_file(file_name)
        parameters = extract_lens_parameters(file_name, lens_type)
        if parameters:
            # Print detailed parameters
            print(f"Model: {file_name}")
            print(f"X Coordinate: {parameters.get('x_coord', '-')}")
            print(f"Y Coordinate: {parameters.get('y_coord', '-')}")
            print(f"Ellipticity: {parameters.get('ellipticity', '-')}")
            print(f"Position Angle: {parameters.get('position_angle', '-')}")
            print(f"Velocity Dispersion: {parameters.get('velocity_dispersion', '-')}")
            print(f"Einstein Radius: {parameters.get('einstein_radius', '-')}")
            print(f"PWI: {parameters.get('pwi', '-')}")
            # Ensure DM Mass is presented in exponential form if it's a float
            dm_mass_formatted = f"{parameters.get('dm_mass', '-')}"
            if isinstance(parameters.get('dm_mass'), float):
                dm_mass_formatted = f"{parameters['dm_mass']:e}"
            print(f"DM Mass: {dm_mass_formatted}")
            print(f"Concentration Parameter: {parameters.get('concentration', '-')}")
            if 'shear_x_coord' in parameters:
                print(f"Shear X Coordinate: {parameters['shear_x_coord']}")
                print(f"Shear Y Coordinate: {parameters['shear_y_coord']}")
                print(f"Shear Position Angle: {parameters['shear_position_angle']}")
                print(f"Shear Strength: {parameters['shear_strength']}")
                print(f"Convergence: {parameters['convergence']}")
            if 'g2_x_coord' in parameters:
                print(f"G2 X Coordinate: {parameters['g2_x_coord']}")
                print(f"G2 Y Coordinate: {parameters['g2_y_coord']}")
                print(f"G2 Ellipticity: {parameters['g2_ellipticity']}")
                print(f"G2 Position Angle: {parameters['g2_position_angle']}")
                print(f"G2 Velocity Dispersion: {parameters['g2_velocity_dispersion']}")
            print()

            # Continue with LaTeX output formatting
            result_text = f"{model_description} & {parameters.get('x_coord', '-'):0.3f} & {parameters.get('y_coord', '-'):0.3f} & {parameters.get('ellipticity', '-'):0.3f} & {parameters.get('position_angle', '-'):0.3f} & {parameters.get('velocity_dispersion', '-')} & {parameters.get('einstein_radius', '-')} & {parameters.get('pwi', '-')} & {dm_mass_formatted} & {parameters.get('concentration', '-')} & - & - \\\\"
            if 'shear_x_coord' in parameters:
                result_text += f"\n - & Shear & - & {parameters['shear_x_coord']:0.3f} & {parameters['shear_y_coord']:0.3f} & - & {parameters['shear_position_angle']:0.3f} & - & - & - & - & - & {parameters['shear_strength']:0.3f} & {parameters['convergence']:0.3f} \\\\"
            if 'g2_x_coord' in parameters:
                result_text += f"\n - & G2 & - & {parameters['g2_x_coord']:0.3f} & {parameters['g2_y_coord']:0.3f} & {parameters['g2_ellipticity']:0.3f} & {parameters['g2_position_angle']:0.3f} & {parameters['g2_velocity_dispersion']} & - & - & - & - & - & - \\\\"
            result_text += "\n\\hline"
            results.append(result_text)
        else:
            result_text = f"{model_description} & No valid data for {lens_type.upper()} \\\\"
            result_text += "\n\\hline"
            results.append(result_text)
    return results

# List of file paths to process
file_list = [
#    "WGDJ0405/SPC/outSP_optresult.dat",
#    "WGDJ0405/SPC/outSPR_optresult.dat",
#    "WGDJ0405/SPFC/outSPF_optresult.dat",
#    "WGDJ0405/SPFC/outSPFR_optresult.dat",
#    "WGDJ0405/PPC/outPP_optresult.dat",
    "WGDJ0405/PPC/outPP_N_optresult.dat",
#    "WGDJ0405/PPC/outPPR_optresult.dat",
    "WGDJ0405/PPC/outPPR_N_optresult.dat",
#    "WGDJ0405/PPFC/outPPF_optresult.dat",
#    "WGDJ0405/PPFC/outPPFR_optresult.dat",
#    "WGDJ0405/NPC/outNP_optresult.dat",
#    "WGDJ0405/NPC/outNPR_optresult.dat",
#    "WGDJ0405/NPFC/outNPF_optresult.dat",
#    "WGDJ0405/NPFC/outNPFR_optresult.dat",
]

# Execute the function and print LaTeX-formatted results
results = process_multiple_files(file_list)
for result in results:
    print(result)


# In[6]:


# Fixed Combined lens parameter table code
# WGD2038

def extract_lens_parameters(filename, lens_type):
    with open(filename, 'r') as file:
        lines = file.readlines()

    lens_key_map = {'s': 'sie', 'p': 'pow', 'n': 'anfw'}
    lens_key = lens_key_map[lens_type.lower()]
    shear_present = 'R' in filename.split('/')[-1]
    g2_present = "G" in filename.split('/')[-1]  # Check if 'G' is in the file name

    last_chi2_index = None
    for index, line in enumerate(lines):
        if "chi^2" in line:
            last_chi2_index = index

    if last_chi2_index is None:
        print(f"Error: 'chi^2' line not found in file: {filename}.")
        return None

    result = {}
    lens_found_count = 0
    first_lens_parameters = {}
    g2_parameters = {}
    found_sie_for_g2 = False

    for line in lines[last_chi2_index:]:
        if "lens   " + lens_key in line:
            parts = line.split()
            try:
                current_params = {
                    'x_coord': float(parts[4]),
                    'y_coord': float(parts[5]),
                    'ellipticity': float(parts[6]),
                    'position_angle': float(parts[7]),
                    'velocity_dispersion': float(parts[3]) if lens_type.lower() == 's' else "-",
                    'einstein_radius': float(parts[8]) if lens_type.lower() == 'p' else "-",
                    'pwi': float(parts[9]) if lens_type.lower() == 'p' else "-",
                    'dm_mass': float(parts[3]) if lens_type.lower() == 'n' else "-",
                    'concentration': float(parts[8]) if lens_type.lower() == 'n' else "-"
                }
                lens_found_count += 1
                if lens_found_count == 1:
                    first_lens_parameters = current_params
            except (ValueError, IndexError):
                print(f"Error: Failed to extract parameters from line: {line} in file: {filename}.")
                return None

        # Check specifically for 'lens   sie' after chi^2 for G2 detection in POW and ANFW models
        if "lens   sie" in line and g2_present:
            parts = line.split()
            try:
                g2_parameters = {
                    'x_coord': float(parts[4]),
                    'y_coord': float(parts[5]),
                    'ellipticity': float(parts[6]),
                    'position_angle': float(parts[7]),
                    'velocity_dispersion': float(parts[3])
                }
                found_sie_for_g2 = True
            except (ValueError, IndexError):
                print(f"Error: Failed to extract G2 parameters from line: {line} in file: {filename}.")
                return None

        if shear_present and "lens   pert" in line:
            parts = line.split()
            try:
                result.update({
                    'shear_x_coord': float(parts[4]),
                    'shear_y_coord': float(parts[5]),
                    'shear_position_angle': float(parts[7]),
                    'shear_strength': float(parts[6]),
                    'convergence': float(parts[9])
                })
            except (ValueError, IndexError):
                print(f"Error: Failed to extract shear parameters from line: {line} in file: {filename}.")
                return None

    result.update(first_lens_parameters)  # Add first lens parameters to result

    # Include G2 parameters in the result if applicable and found
    if g2_present and found_sie_for_g2:
        result.update({
            'g2_x_coord': g2_parameters.get('x_coord'),
            'g2_y_coord': g2_parameters.get('y_coord'),
            'g2_ellipticity': g2_parameters.get('ellipticity'),
            'g2_position_angle': g2_parameters.get('position_angle'),
            'g2_velocity_dispersion': g2_parameters.get('velocity_dispersion')
        })

    if not result:
        print(f"Error: '{lens_key}' line not found after last 'chi^2' in file: {filename}.")
        return None
    return result

def get_lens_type_constraint(file_name):
    """
    Extracts the lens type and constraint from the filename.
    """
    base_name = file_name.split('/')[-1]  # Get the file name part
    lens_type = base_name[3]  # 4th character is lens type
    constraint = base_name[4:6]  # 5th and 6th characters are constraint
    return lens_type.lower(), constraint

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
        lens_type, _ = get_lens_type_constraint(file_name)
        model_description = rename_file(file_name)
        parameters = extract_lens_parameters(file_name, lens_type)
        if parameters:
            # Print detailed parameters
            print(f"Model: {file_name}")
            print(f"X Coordinate: {parameters.get('x_coord', '-')}")
            print(f"Y Coordinate: {parameters.get('y_coord', '-')}")
            print(f"Ellipticity: {parameters.get('ellipticity', '-')}")
            print(f"Position Angle: {parameters.get('position_angle', '-')}")
            print(f"Velocity Dispersion: {parameters.get('velocity_dispersion', '-')}")
            print(f"Einstein Radius: {parameters.get('einstein_radius', '-')}")
            print(f"PWI: {parameters.get('pwi', '-')}")
            # Ensure DM Mass is presented in exponential form if it's a float
            dm_mass_formatted = f"{parameters.get('dm_mass', '-')}"
            if isinstance(parameters.get('dm_mass'), float):
                dm_mass_formatted = f"{parameters['dm_mass']:e}"
            print(f"DM Mass: {dm_mass_formatted}")
            print(f"Concentration Parameter: {parameters.get('concentration', '-')}")
            if 'shear_x_coord' in parameters:
                print(f"Shear X Coordinate: {parameters['shear_x_coord']}")
                print(f"Shear Y Coordinate: {parameters['shear_y_coord']}")
                print(f"Shear Position Angle: {parameters['shear_position_angle']}")
                print(f"Shear Strength: {parameters['shear_strength']}")
                print(f"Convergence: {parameters['convergence']}")
            if 'g2_x_coord' in parameters:
                print(f"G2 X Coordinate: {parameters['g2_x_coord']}")
                print(f"G2 Y Coordinate: {parameters['g2_y_coord']}")
                print(f"G2 Ellipticity: {parameters['g2_ellipticity']}")
                print(f"G2 Position Angle: {parameters['g2_position_angle']}")
                print(f"G2 Velocity Dispersion: {parameters['g2_velocity_dispersion']}")
            print()

            # Continue with LaTeX output formatting
            result_text = f"{model_description} & {parameters.get('x_coord', '-'):0.3f} & {parameters.get('y_coord', '-'):0.3f} & {parameters.get('ellipticity', '-'):0.3f} & {parameters.get('position_angle', '-'):0.3f} & {parameters.get('velocity_dispersion', '-')} & {parameters.get('einstein_radius', '-')} & {parameters.get('pwi', '-')} & {dm_mass_formatted} & {parameters.get('concentration', '-')} & - & - \\\\"
            if 'shear_x_coord' in parameters:
                result_text += f"\n - & Shear & - & {parameters['shear_x_coord']:0.3f} & {parameters['shear_y_coord']:0.3f} & - & {parameters['shear_position_angle']:0.3f} & - & - & - & - & - & {parameters['shear_strength']:0.3f} & {parameters['convergence']:0.3f} \\\\"
            if 'g2_x_coord' in parameters:
                result_text += f"\n - & G2 & - & {parameters['g2_x_coord']:0.3f} & {parameters['g2_y_coord']:0.3f} & {parameters['g2_ellipticity']:0.3f} & {parameters['g2_position_angle']:0.3f} & {parameters['g2_velocity_dispersion']} & - & - & - & - & - & - \\\\"
            result_text += "\n\\hline"
            results.append(result_text)
        else:
            result_text = f"{model_description} & No valid data for {lens_type.upper()} \\\\"
            result_text += "\n\\hline"
            results.append(result_text)
    return results

# List of file paths to process
file_list = [
#    "WGD2038/SPC/outSP_optresult.dat",
#    "WGD2038/SPC/outSPR_optresult.dat",
#    "WGD2038/SPFC/outSPF_optresult.dat",
#    "WGD2038/SPFC/outSPFR_optresult.dat",
#    "WGD2038/PPC/outPP_optresult.dat",
    "WGD2038/PPC/outPP_N_Free_PWI_optresult.dat",
    "WGD2038/PPC/outPPR_N_Free_PWI_optresult.dat",
#    "WGD2038/PPC/outPP_N_Free_PWI_Free_PA_optresult.dat",
#    "WGD2038/PPC/outPPR_N_Free_PWI_Free_PA_optresult.dat",
#    "WGD2038/PPC/outPPR_optresult.dat",
#    "WGD2038/PPFC/outPPF_optresult.dat",
#    "WGD2038/PPFC/outPPFR_optresult.dat",
#    "WGD2038/NPC/outNP_optresult.dat",
#    "WGD2038/NPC/outNPR_optresult.dat",
#    "WGD2038/NPFC/outNPF_optresult.dat",
#    "WGD2038/NPFC/outNPFR_optresult.dat",
]

# Execute the function and print LaTeX-formatted results
results = process_multiple_files(file_list)
for result in results:
    print(result)

