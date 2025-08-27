#!/usr/bin/env python
# coding: utf-8

# In[3]:


def extract_sie_parameters(filename):
    """
    Extract the SHEAR lens parameters from the last "lens pert" line that appears in the file.
    """
    with open(filename, 'r') as file:
        lines = file.readlines()

    last_sie_index = None
    for index, line in enumerate(lines):
        if "lens   pert" in line:
            last_sie_index = index

    if last_sie_index is None:
        print(f"Error: 'lens pert' line not found in file: {filename}.")
        return None

    line = lines[last_sie_index]
    parts = line.split()
    try:
        lens_redshift = float(parts[3])
        source_redshift = float(parts[2])
        x_coord = float(parts[4])
        y_coord = float(parts[5])
        shear_strength = float(parts[6])
        shear_angle = float(parts[7])
        convergence = float(parts[9])
        return (lens_redshift, source_redshift, x_coord, y_coord, shear_strength, shear_angle, convergence)
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
            source_redshift, lens_redshift, x_coord, y_coord, shear_strength, shear_angle, convergence = sie_parameters
            result_text = f"{model_description} & {lens_redshift:.3f} & {source_redshift:.3f} & {x_coord:.3f} & {y_coord:.3f} & {shear_strength:.3f} & {shear_angle:.3f} & {convergence:.3f} \\\\"
            results.append(result_text)
            print(f"Model: {file_name}")
            print(f"Lens Redshift: {lens_redshift:.3f}")
            print(f"Source Redshift: {source_redshift:.3f}")
            print(f"X Coordinate: {x_coord:.3f}")
            print(f"Y Coordinate: {y_coord:.3f}")
            print(f"Shear Strength: {shear_strength:.3f}")
            print(f"Shear Position Angle: {shear_angle:.3f}")
            print(f"Convergence: {convergence:.3f}\n")
        else:
            results.append(f"{model_description} Failed to extract SIE parameters.")
            print(f"Model: {file_name} Failed to extract SIE parameters.\n")
    return results

# List of file paths to process
file_list = [
    "RXJ0911/SPC/outSPR_optresult.dat",
    "RXJ0911/SPC/outSPGR_optresult.dat",
    "RXJ0911/SPFC/outSPFR_optresult.dat",
    "RXJ0911/SPFC/outSPFGR_optresult.dat",
    "RXJ0911/PPC/outPPR_optresult.dat",
    "RXJ0911/PPC/outPPGR_optresult.dat",
    "RXJ0911/PPFC/outPPFR_optresult.dat",
    "RXJ0911/PPFC/outPPFGR_optresult.dat",
    "RXJ0911/NPC/outNPR_optresult.dat",
    "RXJ0911/NPC/outNPGR_optresult.dat",
    "RXJ0911/NPFC/outNPFR_optresult.dat",
    "RXJ0911/NPFC/outNPFGR_optresult.dat",
    "PSJ1606/SPC/outSPR_optresult.dat",
    "PSJ1606/SPC/outSPGR_optresult.dat",
    "PSJ1606/SPFC/outSPFR_optresult.dat",
    "PSJ1606/SPFC/outSPFGR_optresult.dat",
    "PSJ1606/PPC/outPPR_optresult.dat",
    "PSJ1606/PPC/outPPGR_optresult.dat",
    "PSJ1606/PPFC/outPPFR_optresult.dat",
    "PSJ1606/PPFC/outPPFGR_optresult.dat",
    "PSJ1606/NPC/outNPR_optresult.dat",
    "PSJ1606/NPC/outNPGR_optresult.dat",
    "PSJ1606/NPFC/outNPFR_optresult.dat",
    "PSJ1606/NPFC/outNPFGR_optresult.dat",
    "WFI2033/SPR/outSPR_optresult.dat",
    "WFI2033/SPGR/outSPGR_optresult.dat",
    "WFI2033/SPFR/outSPFR_optresult.dat",
    "WFI2033/SPFGR/outSPFGR_optresult.dat",
    "WFI2033/PPR/outPPR_optresult.dat",
    "WFI2033/PPGR/outPPGR_optresult.dat",
    "WFI2033/PPFR/outPPFR_optresult.dat",
    "WFI2033/PPFGR/outPPFGR_optresult.dat",
    "WFI2033/NPR/outNPR_optresult.dat",
    "WFI2033/NPGR/outNPGR_optresult.dat",
    "WFI2033/NPFR/outNPFR_optresult.dat",
    "WFI2033/NPFGR/outNPFGR_optresult.dat",
    "SDSSJ1330/SPC/outSPR_optresult.dat",
    "SDSSJ1330/SPFC/outSPFR_optresult.dat",
    "SDSSJ1330/PPC/outPPR_optresult.dat",
    "SDSSJ1330/PPFC/outPPFR_optresult.dat",
    "SDSSJ1330/NPC/outNPR_optresult.dat",
    "SDSSJ1330/NPFC/outNPFR_optresult.dat",
    "WFI2026/SPC/outSPR_optresult.dat",
    "WFI2026/SPFC/outSPFR_optresult.dat",
    "WFI2026/PPC/outPPR_optresult.dat",
    "WFI2026/PPFC/outPPFR_optresult.dat",
    "WFI2026/NPC/outNPR_optresult.dat",
    "WFI2026/NPFC/outNPFR_optresult.dat",
    "WGDJ0405/SPC/outSPR_optresult.dat",
    "WGDJ0405/SPFC/outSPFR_optresult.dat",
    "WGDJ0405/PPC/outPPR_optresult.dat",
    "WGDJ0405/PPFC/outPPFR_optresult.dat",
    "WGDJ0405/NPC/outNPR_optresult.dat",
    "WGDJ0405/NPFC/outNPFR_optresult.dat",
    "WGD2038/SPC/outSPR_optresult.dat",
    "WGD2038/SPFC/outSPFR_optresult.dat",
    "WGD2038/PPC/outPPR_optresult.dat",
    "WGD2038/PPFC/outPPFR_optresult.dat",
    "WGD2038/NPC/outNPR_optresult.dat",
    "WGD2038/NPFC/outNPFR_optresult.dat",
]

# Execute the function and print LaTeX-formatted results
latex_results = process_multiple_files(file_list)
for result in latex_results:
    print(result)


# In[ ]:




