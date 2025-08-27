#!/usr/bin/env python
# coding: utf-8

# In[3]:


def extract_sie_parameters(filename):
    """
    Extract the NFW lens parameters from the last "lens anfw" line that appears in the file.
    """
    with open(filename, 'r') as file:
        lines = file.readlines()

    last_sie_index = None
    for index, line in enumerate(lines):
        if "lens   anfw" in line:
            last_sie_index = index

    if last_sie_index is None:
        print(f"Error: 'lens anfw' line not found in file: {filename}.")
        return None

    line = lines[last_sie_index]
    parts = line.split()
    try:
        lens_redshift = float(parts[2])
        dm_mass = float(parts[3])
        x_coord = float(parts[4])
        y_coord = float(parts[5])
        ellipticity = float(parts[6])
        position_angle = float(parts[7])
        concentration_para = float(parts[8])
        return (lens_redshift, dm_mass, x_coord, y_coord, ellipticity, position_angle, concentration_para)
    except (ValueError, IndexError):
        print(f"Error: Failed to extract numbers from line: {line} in file: {filename}.")
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
        constraint = "QSO Pos + FR"
    elif key[1] == "P":
        constraint = "QSO Pos"

    extras = []
    if "G" in key:
        extras.append("G2")
    if "R" in key:
        extras.append("Shear")

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
            lens_redshift, dm_mass, x_coord, y_coord, ellipticity, position_angle, concentration_para = sie_parameters
            result_text = f"{model_description} & {lens_redshift:.3f} & {dm_mass:e} & {x_coord:.3f} & {y_coord:.3f} & {ellipticity:.3f} & {position_angle:.3f} & {concentration_para:.3f} \\\\"
            results.append(result_text)
            print(f"Model: {file_name}")
            print(f"Lens Redshift: {lens_redshift:.3f}")
            print(f"DM Mass: {dm_mass:e}")
            print(f"X Coordinate: {x_coord:.3f}")
            print(f"Y Coordinate: {y_coord:.3f}")
            print(f"Ellipticity: {ellipticity:.3f}")
            print(f"Position Angle: {position_angle:.3f}")
            print(f"Concentration Parameter: {concentration_para:.3f}\n")
        else:
            results.append(f"{model_description} Failed to extract SIE parameters.")
            print(f"Model: {file_name} Failed to extract SIE parameters.\n")
    return results

# List of file paths to process
file_list = [
    "RXJ0911/NPC/outNP_optresult.dat",
    "RXJ0911/NPC/outNPR_optresult.dat",
    "RXJ0911/NPC/outNPG_optresult.dat",
    "RXJ0911/NPC/outNPGR_optresult.dat",
    "RXJ0911/NPFC/outNPF_optresult.dat",
    "RXJ0911/NPFC/outNPFR_optresult.dat",
    "RXJ0911/NPFC/outNPFG_optresult.dat",
    "RXJ0911/NPFC/outNPFGR_optresult.dat",
    "PSJ1606/NPC/outNP_optresult.dat",
    "PSJ1606/NPC/outNPR_optresult.dat",
    "PSJ1606/NPC/outNPG_optresult.dat",
    "PSJ1606/NPC/outNPGR_optresult.dat",
    "PSJ1606/NPFC/outNPF_optresult.dat",
    "PSJ1606/NPFC/outNPFR_optresult.dat",
    "PSJ1606/NPFC/outNPFG_optresult.dat",
    "PSJ1606/NPFC/outNPFGR_optresult.dat",
    "WFI2033/NP/outNP_optresult.dat",
    "WFI2033/NPR/outNPR_optresult.dat",
    "WFI2033/NPG/outNPG_optresult.dat",
    "WFI2033/NPGR/outNPGR_optresult.dat",
    "WFI2033/NPF/outNPF_optresult.dat",
    "WFI2033/NPFR/outNPFR_optresult.dat",
    "WFI2033/NPFG/outNPFG_optresult.dat",
    "WFI2033/NPFGR/outNPFGR_optresult.dat",
    "SDSSJ1330/NPC/outNP_optresult.dat",
    "SDSSJ1330/NPC/outNPR_optresult.dat",
    "SDSSJ1330/NPFC/outNPF_optresult.dat",
    "SDSSJ1330/NPFC/outNPFR_optresult.dat",
    "WFI2026/NPC/outNP_optresult.dat",
    "WFI2026/NPC/outNPR_optresult.dat",
    "WFI2026/NPFC/outNPF_optresult.dat",
    "WFI2026/NPFC/outNPFR_optresult.dat",
    "WGDJ0405/NPC/outNP_optresult.dat",
    "WGDJ0405/NPC/outNPR_optresult.dat",
    "WGDJ0405/NPFC/outNPF_optresult.dat",
    "WGDJ0405/NPFC/outNPFR_optresult.dat",
    "WGD2038/NPC/outNP_optresult.dat",
    "WGD2038/NPC/outNPR_optresult.dat",
    "WGD2038/NPFC/outNPF_optresult.dat",
    "WGD2038/NPFC/outNPFR_optresult.dat"
]

# Execute the function and print LaTeX-formatted results
latex_results = process_multiple_files(file_list)
for result in latex_results:
    print(result)


# In[ ]:




