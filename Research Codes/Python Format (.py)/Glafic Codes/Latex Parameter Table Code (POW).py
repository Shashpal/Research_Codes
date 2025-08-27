#!/usr/bin/env python
# coding: utf-8

# In[1]:


def extract_sie_parameters(filename):
    """
    Extract the POW lens parameters from the last "lens pow" line that appears in the file.
    """
    with open(filename, 'r') as file:
        lines = file.readlines()

    last_sie_index = None
    for index, line in enumerate(lines):
        if "lens   pow" in line:
            last_sie_index = index

    if last_sie_index is None:
        print(f"Error: 'lens pow' line not found in file: {filename}.")
        return None

    line = lines[last_sie_index]
    parts = line.split()
    try:
        lens_redshift = float(parts[3])
        source_redshift = float(parts[2])
        x_coord = float(parts[4])
        y_coord = float(parts[5])
        ellipticity = float(parts[6])
        position_angle = float(parts[7])
        einstein_radius = float(parts[8])
        pwi = float(parts[9])
        return (lens_redshift, source_redshift, x_coord, y_coord, ellipticity, position_angle, einstein_radius, pwi)
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
            source_redshift, lens_redshift, x_coord, y_coord, ellipticity, position_angle, einstein_radius, pwi = sie_parameters
            result_text = f"{model_description} & {lens_redshift:.3f} & {source_redshift:.3f} & {x_coord:.3f} & {y_coord:.3f} & {ellipticity:.3f} & {position_angle:.3f} & {einstein_radius:.3f} & {pwi:.3f} \\\\"
            results.append(result_text)
            print(f"Model: {file_name}")
            print(f"Lens Redshift: {lens_redshift:.3f}")
            print(f"Source Redshift: {source_redshift:.3f}")
            print(f"X Coordinate: {x_coord:.3f}")
            print(f"Y Coordinate: {y_coord:.3f}")
            print(f"Ellipticity: {ellipticity:.3f}")
            print(f"Position Angle: {position_angle:.3f}")
            print(f"Einstein Radius: {einstein_radius:.3f}")
            print(f"PWI: {pwi:.3f}\n")
        else:
            results.append(f"{model_description} Failed to extract SIE parameters.")
            print(f"Model: {file_name} Failed to extract SIE parameters.\n")
    return results

# List of file paths to process
file_list = [
    "RXJ0911/PPC/outPP_optresult.dat",
    "RXJ0911/PPC/outPPR_optresult.dat",
    "RXJ0911/PPC/outPPG_optresult.dat",
    "RXJ0911/PPC/outPPGR_optresult.dat",
    "RXJ0911/PPFC/outPPF_optresult.dat",
    "RXJ0911/PPFC/outPPFR_optresult.dat",
    "RXJ0911/PPFC/outPPFG_optresult.dat",
    "RXJ0911/PPFC/outPPFGR_optresult.dat",
    "PSJ1606/PPC/outPP_optresult.dat",
    "PSJ1606/PPC/outPPR_optresult.dat",
    "PSJ1606/PPC/outPPG_optresult.dat",
    "PSJ1606/PPC/outPPGR_optresult.dat",
    "PSJ1606/PPFC/outPPF_optresult.dat",
    "PSJ1606/PPFC/outPPFR_optresult.dat",
    "PSJ1606/PPFC/outPPFG_optresult.dat",
    "PSJ1606/PPFC/outPPFGR_optresult.dat",
    "WFI2033/PP/outPP_optresult.dat",
    "WFI2033/PPR/outPPR_optresult.dat",
    "WFI2033/PPG/outPPG_optresult.dat",
    "WFI2033/PPGR/outPPGR_optresult.dat",
    "WFI2033/PPF/outPPF_optresult.dat",
    "WFI2033/PPFR/outPPFR_optresult.dat",
    "WFI2033/PPFG/outPPFG_optresult.dat",
    "WFI2033/PPFGR/outPPFGR_optresult.dat",
    "SDSSJ1330/PPC/outPP_optresult.dat",
    "SDSSJ1330/PPC/outPPR_optresult.dat",
    "SDSSJ1330/PPFC/outPPF_optresult.dat",
    "SDSSJ1330/PPFC/outPPFR_optresult.dat",
    "WFI2026/PPC/outPP_optresult.dat",
    "WFI2026/PPC/outPPR_optresult.dat",
    "WFI2026/PPFC/outPPF_optresult.dat",
    "WFI2026/PPFC/outPPFR_optresult.dat",
    "WGDJ0405/PPC/outPP_optresult.dat",
    "WGDJ0405/PPC/outPPR_optresult.dat",
    "WGDJ0405/PPFC/outPPF_optresult.dat",
    "WGDJ0405/PPFC/outPPFR_optresult.dat",
    "WGD2038/PPC/outPP_optresult.dat",
    "WGD2038/PPC/outPPR_optresult.dat",
    "WGD2038/PPFC/outPPF_optresult.dat",
    "WGD2038/PPFC/outPPFR_optresult.dat",
]

# Execute the function and print LaTeX-formatted results
latex_results = process_multiple_files(file_list)
for result in latex_results:
    print(result)


# In[ ]:




