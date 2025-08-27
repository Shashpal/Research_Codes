#!/usr/bin/env python
# coding: utf-8

# In[2]:


# RXJ0911 lens centre offset code
import math
from astropy.cosmology import Planck18 as cosmo

def extract_coordinates(filename):
    """
    Extract the (x, y) coordinates for the first and second mentioned profiles after the latest "chi^2" line.
    The first profile is labeled as G, and the second profile is labeled as G2.
    """
    x_coordinate_g = None
    y_coordinate_g = None
    x_coordinate_g2 = None
    y_coordinate_g2 = None

    try:
        with open(filename, 'r') as file:
            lines = file.readlines()

        # Find the index of the latest "chi^2" line
        latest_chi2_index = None
        for i, line in enumerate(lines):
            if "chi^2" in line:
                latest_chi2_index = i

        if latest_chi2_index is None:
            return None, None, None, None

        # Process lines after the latest "chi^2" line
        for line in lines[latest_chi2_index + 1:]:
            if "lens   sie" in line or "lens   pow" in line or "lens   anfw" in line:
                # Extract numerical values using a fixed format approach
                parts = [value for value in line.split() if value]

                try:
                    # Extract the x and y coordinates (3rd and 4th parameters)
                    x_coord = float(parts[4])  # 3rd parameter
                    y_coord = float(parts[5])  # 4th parameter

                    # Assign G and G2 based on the order of appearance
                    if x_coordinate_g is None and y_coordinate_g is None:
                        x_coordinate_g = x_coord
                        y_coordinate_g = y_coord
                    elif x_coordinate_g2 is None and y_coordinate_g2 is None:
                        x_coordinate_g2 = x_coord
                        y_coordinate_g2 = y_coord
                        break  # Stop after finding G and G2
                except (ValueError, IndexError):
                    return x_coordinate_g, y_coordinate_g, None, None

    except FileNotFoundError:
        return None, None, None, None
    except Exception:
        return None, None, None, None

    return x_coordinate_g, y_coordinate_g, x_coordinate_g2, y_coordinate_g2


def rename_file(file_path):
    """
    Rename a file based on its naming conventions.
    """
    try:
        parts = file_path.split("/")
        lens_system = parts[0]
        file_name = parts[-1]
        key = file_name.split("_optresult")[0][3:]
        profile_map = {"S": "SIE", "P": "POW", "N": "NFW"}
        profile = profile_map.get(key[0], "Unknown")
        return f"{lens_system} {profile} {key}"
    except IndexError:
        return "Unknown Description"


def calculate_offset(x, y, observed_x, observed_y):
    """
    Calculate the offset in absolute scale.
    """
    return abs(x - observed_x), abs(y - observed_y)


def calculate_total_magnitude(offset_x, offset_y):
    """
    Calculate the total magnitude of the offset using the formula:
    sqrt((x1 - x2)**2 + (y1 - y2)**2).
    """
    return math.sqrt(offset_x**2 + offset_y**2)


def calculate_physical_offset(magnitude_arcsec, redshift):
    """
    Calculate the physical offset in meters using the formula L = r * θ.
    θ is the magnitude in radians, and r is the angular diameter distance.
    """
    theta_radians = math.radians(magnitude_arcsec / 3600.0)
    r_meters = cosmo.angular_diameter_distance(redshift).to('m').value
    l_meters = r_meters * theta_radians
    l_pc = l_meters / 3.086e16
    return l_pc


def calculate_reference_scale(redshift):
    """
    Calculate the arcsecond offset for a physical distance of 10 pc.
    """
    distance_10pc_meters = 10 * 3.086e16
    r_meters = cosmo.angular_diameter_distance(redshift).to('m').value
    theta_radians = distance_10pc_meters / r_meters
    theta_arcseconds = math.degrees(theta_radians) * 3600.0
    return theta_arcseconds


def process_multiple_files(file_list, output_filename):
    """
    Process a list of files, rename them, and write the extracted results including offsets to a file,
    while also printing the results to the console.
    """
    observed_g = (0, 0)
    observed_g2 = (-0.767, 0.657)
    redshift = 0.769
    reference_scale = calculate_reference_scale(redshift)

    with open(output_filename, 'w') as output_file:
        # Write reference scale at the top
        header = f"Reference Scale: (x, y) = ({reference_scale:.4f}, {reference_scale:.4f}) arcsec for 10 pc\n\n"
        print(header.strip())
        output_file.write(header)

        for file_name in file_list:
            # Rename the file based on naming conventions
            description = rename_file(file_name)

            # Extract x and y coordinates for G and G2
            x_g, y_g, x_g2, y_g2 = extract_coordinates(file_name)

            if x_g is not None and y_g is not None:
                offset_x_g, offset_y_g = calculate_offset(x_g, y_g, observed_g[0], observed_g[1])
                total_magnitude_g = calculate_total_magnitude(offset_x_g, offset_y_g)
                physical_offset_g = calculate_physical_offset(total_magnitude_g, redshift)
                result_g = (
                    f"{description} : G (x, y) = ({x_g:.4f}, {y_g:.4f}), "
                    f"Offset = ({offset_x_g:.4f}, {offset_y_g:.4f}), "
                    f"Magnitude = {total_magnitude_g:.3f}, Physical Offset = {physical_offset_g:.0f} pc"
                )
                print(result_g)
                output_file.write(result_g + "\n")
            else:
                result_g = f"{description}: Failed to extract G (x, y)."
                print(result_g)
                output_file.write(result_g + "\n")

            if x_g2 is not None and y_g2 is not None:
                offset_x_g2, offset_y_g2 = calculate_offset(x_g2, y_g2, observed_g2[0], observed_g2[1])
                total_magnitude_g2 = calculate_total_magnitude(offset_x_g2, offset_y_g2)
                physical_offset_g2 = calculate_physical_offset(total_magnitude_g2, redshift)
                result_g2 = (
                    f"{description} : G2 (x, y) = ({x_g2:.4f}, {y_g2:.4f}), "
                    f"Offset = ({offset_x_g2:.4f}, {offset_y_g2:.4f}), "
                    f"Magnitude = {total_magnitude_g2:.3f}, Physical Offset = {physical_offset_g2:.0f} pc"
                )
                print(result_g2)
                output_file.write(result_g2 + "\n")


# List of file paths to process
file_list = [
    "RXJ0911/SPC/outSP_optresult.dat",
    "RXJ0911/SPC/outSPG_optresult.dat",
    "RXJ0911/SPC/outSPR_optresult.dat",
    "RXJ0911/SPC/outSPGR_optresult.dat",
    "RXJ0911/SPFC/outSPF_optresult.dat",
    "RXJ0911/SPFC/outSPFG_optresult.dat",
    "RXJ0911/SPFC/outSPFR_optresult.dat",
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

# Process the files and write results to a text file
process_multiple_files(file_list, "RXJ0911_lens_centre_offset_results_gau_priored.txt")


# In[3]:


# PSJ1606 lens centre offset code
import math
from astropy.cosmology import Planck18 as cosmo

def extract_coordinates(filename):
    """
    Extract the (x, y) coordinates for the first and second mentioned profiles after the latest "chi^2" line.
    The first profile is labeled as G, and the second profile is labeled as G2.
    """
    x_coordinate_g = None
    y_coordinate_g = None
    x_coordinate_g2 = None
    y_coordinate_g2 = None

    try:
        with open(filename, 'r') as file:
            lines = file.readlines()

        # Find the index of the latest "chi^2" line
        latest_chi2_index = None
        for i, line in enumerate(lines):
            if "chi^2" in line:
                latest_chi2_index = i

        if latest_chi2_index is None:
            return None, None, None, None

        # Process lines after the latest "chi^2" line
        for line in lines[latest_chi2_index + 1:]:
            if "lens   sie" in line or "lens   pow" in line or "lens   anfw" in line:
                # Extract numerical values using a fixed format approach
                parts = [value for value in line.split() if value]

                try:
                    # Extract the x and y coordinates (3rd and 4th parameters)
                    x_coord = float(parts[4])  # 3rd parameter
                    y_coord = float(parts[5])  # 4th parameter

                    # Assign G and G2 based on the order of appearance
                    if x_coordinate_g is None and y_coordinate_g is None:
                        x_coordinate_g = x_coord
                        y_coordinate_g = y_coord
                    elif x_coordinate_g2 is None and y_coordinate_g2 is None:
                        x_coordinate_g2 = x_coord
                        y_coordinate_g2 = y_coord
                        break  # Stop after finding G and G2
                except (ValueError, IndexError):
                    return x_coordinate_g, y_coordinate_g, None, None

    except FileNotFoundError:
        return None, None, None, None
    except Exception:
        return None, None, None, None

    return x_coordinate_g, y_coordinate_g, x_coordinate_g2, y_coordinate_g2


def rename_file(file_path):
    """
    Rename a file based on its naming conventions.
    """
    try:
        parts = file_path.split("/")
        lens_system = parts[0]
        file_name = parts[-1]
        key = file_name.split("_optresult")[0][3:]
        profile_map = {"S": "SIE", "P": "POW", "N": "NFW"}
        profile = profile_map.get(key[0], "Unknown")
        return f"{lens_system} {profile} {key}"
    except IndexError:
        return "Unknown Description"


def calculate_offset(x, y, observed_x, observed_y):
    """
    Calculate the offset in absolute scale.
    """
    return abs(x - observed_x), abs(y - observed_y)


def calculate_total_magnitude(offset_x, offset_y):
    """
    Calculate the total magnitude of the offset using the formula:
    sqrt((x1 - x2)**2 + (y1 - y2)**2).
    """
    return math.sqrt(offset_x**2 + offset_y**2)


def calculate_physical_offset(magnitude_arcsec, redshift):
    """
    Calculate the physical offset in meters using the formula L = r * θ.
    θ is the magnitude in radians, and r is the angular diameter distance.
    """
    theta_radians = math.radians(magnitude_arcsec / 3600.0)
    r_meters = cosmo.angular_diameter_distance(redshift).to('m').value
    l_meters = r_meters * theta_radians
    l_pc = l_meters / 3.086e16
    return l_pc


def calculate_reference_scale(redshift):
    """
    Calculate the arcsecond offset for a physical distance of 10 pc.
    """
    distance_10pc_meters = 10 * 3.086e16
    r_meters = cosmo.angular_diameter_distance(redshift).to('m').value
    theta_radians = distance_10pc_meters / r_meters
    theta_arcseconds = math.degrees(theta_radians) * 3600.0
    return theta_arcseconds


def process_multiple_files(file_list, output_filename):
    """
    Process a list of files, rename them, and write the extracted results including offsets to a file,
    while also printing the results to the console.
    """
    observed_g = (0, 0)
    observed_g2 = (-0.307, -1.153)
    redshift = 0.3
    reference_scale = calculate_reference_scale(redshift)

    with open(output_filename, 'w') as output_file:
        # Write reference scale at the top
        header = f"Reference Scale: (x, y) = ({reference_scale:.4f}, {reference_scale:.4f}) arcsec for 10 pc\n\n"
        print(header.strip())
        output_file.write(header)

        for file_name in file_list:
            # Rename the file based on naming conventions
            description = rename_file(file_name)

            # Extract x and y coordinates for G and G2
            x_g, y_g, x_g2, y_g2 = extract_coordinates(file_name)

            if x_g is not None and y_g is not None:
                offset_x_g, offset_y_g = calculate_offset(x_g, y_g, observed_g[0], observed_g[1])
                total_magnitude_g = calculate_total_magnitude(offset_x_g, offset_y_g)
                physical_offset_g = calculate_physical_offset(total_magnitude_g, redshift)
                result_g = (
                    f"{description} : G (x, y) = ({x_g:.4f}, {y_g:.4f}), "
                    f"Offset = ({offset_x_g:.4f}, {offset_y_g:.4f}), "
                    f"Magnitude = {total_magnitude_g:.3f}, Physical Offset = {physical_offset_g:.0f} pc"
                )
                print(result_g)
                output_file.write(result_g + "\n")
            else:
                result_g = f"{description}: Failed to extract G (x, y)."
                print(result_g)
                output_file.write(result_g + "\n")

            if x_g2 is not None and y_g2 is not None:
                offset_x_g2, offset_y_g2 = calculate_offset(x_g2, y_g2, observed_g2[0], observed_g2[1])
                total_magnitude_g2 = calculate_total_magnitude(offset_x_g2, offset_y_g2)
                physical_offset_g2 = calculate_physical_offset(total_magnitude_g2, redshift)
                result_g2 = (
                    f"{description} : G2 (x, y) = ({x_g2:.4f}, {y_g2:.4f}), "
                    f"Offset = ({offset_x_g2:.4f}, {offset_y_g2:.4f}), "
                    f"Magnitude = {total_magnitude_g2:.3f}, Physical Offset = {physical_offset_g2:.0f} pc"
                )
                print(result_g2)
                output_file.write(result_g2 + "\n")


# List of file paths to process
file_list = [
    "PSJ1606/SPC/outSP_optresult.dat",
    "PSJ1606/SPC/outSPG_optresult.dat",
    "PSJ1606/SPC/outSPR_optresult.dat",
    "PSJ1606/SPC/outSPGR_optresult.dat",
    "PSJ1606/SPFC/outSPF_optresult.dat",
    "PSJ1606/SPFC/outSPFG_optresult.dat",
    "PSJ1606/SPFC/outSPFR_optresult.dat",
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

# Process the files and write results to a text file
process_multiple_files(file_list, "PSJ1606_lens_centre_offset_results_gau_priored.txt")


# In[1]:


# WFI2033 lens centre offset code
import math
from astropy.cosmology import Planck18 as cosmo

def extract_coordinates(filename):
    """
    Extract the (x, y) coordinates for the first and second mentioned profiles after the latest "chi^2" line.
    The first profile is labeled as G, and the second profile is labeled as G2.
    """
    x_coordinate_g = None
    y_coordinate_g = None
    x_coordinate_g2 = None
    y_coordinate_g2 = None

    try:
        with open(filename, 'r') as file:
            lines = file.readlines()

        # Find the index of the latest "chi^2" line
        latest_chi2_index = None
        for i, line in enumerate(lines):
            if "chi^2" in line:
                latest_chi2_index = i

        if latest_chi2_index is None:
            return None, None, None, None

        # Process lines after the latest "chi^2" line
        for line in lines[latest_chi2_index + 1:]:
            if "lens   sie" in line or "lens   pow" in line or "lens   anfw" in line:
                # Extract numerical values using a fixed format approach
                parts = [value for value in line.split() if value]

                try:
                    # Extract the x and y coordinates (3rd and 4th parameters)
                    x_coord = float(parts[4])  # 3rd parameter
                    y_coord = float(parts[5])  # 4th parameter

                    # Assign G and G2 based on the order of appearance
                    if x_coordinate_g is None and y_coordinate_g is None:
                        x_coordinate_g = x_coord
                        y_coordinate_g = y_coord
                    elif x_coordinate_g2 is None and y_coordinate_g2 is None:
                        x_coordinate_g2 = x_coord
                        y_coordinate_g2 = y_coord
                        break  # Stop after finding G and G2
                except (ValueError, IndexError):
                    return x_coordinate_g, y_coordinate_g, None, None

    except FileNotFoundError:
        return None, None, None, None
    except Exception:
        return None, None, None, None

    return x_coordinate_g, y_coordinate_g, x_coordinate_g2, y_coordinate_g2


def rename_file(file_path):
    """
    Rename a file based on its naming conventions.
    """
    try:
        parts = file_path.split("/")
        lens_system = parts[0]
        file_name = parts[-1]
        key = file_name.split("_optresult")[0][3:]
        profile_map = {"S": "SIE", "P": "POW", "N": "NFW"}
        profile = profile_map.get(key[0], "Unknown")
        return f"{lens_system} {profile} {key}"
    except IndexError:
        return "Unknown Description"


def calculate_offset(x, y, observed_x, observed_y):
    """
    Calculate the offset in absolute scale.
    """
    return abs(x - observed_x), abs(y - observed_y)


def calculate_total_magnitude(offset_x, offset_y):
    """
    Calculate the total magnitude of the offset using the formula:
    sqrt((x1 - x2)**2 + (y1 - y2)**2).
    """
    return math.sqrt(offset_x**2 + offset_y**2)


def calculate_physical_offset(magnitude_arcsec, redshift):
    """
    Calculate the physical offset in meters using the formula L = r * θ.
    θ is the magnitude in radians, and r is the angular diameter distance.
    """
    theta_radians = math.radians(magnitude_arcsec / 3600.0)
    r_meters = cosmo.angular_diameter_distance(redshift).to('m').value
    l_meters = r_meters * theta_radians
    l_pc = l_meters / 3.086e16
    return l_pc


def calculate_reference_scale(redshift):
    """
    Calculate the arcsecond offset for a physical distance of 10 pc.
    """
    distance_10pc_meters = 10 * 3.086e16
    r_meters = cosmo.angular_diameter_distance(redshift).to('m').value
    theta_radians = distance_10pc_meters / r_meters
    theta_arcseconds = math.degrees(theta_radians) * 3600.0
    return theta_arcseconds


def process_multiple_files(file_list, output_filename):
    """
    Process a list of files, rename them, and write the extracted results including offsets to a file,
    while also printing the results to the console.
    """
    observed_g = (0, 0)
    observed_g2 = (0.245, 2.037)
    redshift = 0.661
    reference_scale = calculate_reference_scale(redshift)

    with open(output_filename, 'w') as output_file:
        # Write reference scale at the top
        header = f"Reference Scale: (x, y) = ({reference_scale:.4f}, {reference_scale:.4f}) arcsec for 10 pc\n\n"
        print(header.strip())
        output_file.write(header)

        for file_name in file_list:
            # Rename the file based on naming conventions
            description = rename_file(file_name)

            # Extract x and y coordinates for G and G2
            x_g, y_g, x_g2, y_g2 = extract_coordinates(file_name)

            if x_g is not None and y_g is not None:
                offset_x_g, offset_y_g = calculate_offset(x_g, y_g, observed_g[0], observed_g[1])
                total_magnitude_g = calculate_total_magnitude(offset_x_g, offset_y_g)
                physical_offset_g = calculate_physical_offset(total_magnitude_g, redshift)
                result_g = (
                    f"{description} : G (x, y) = ({x_g:.4f}, {y_g:.4f}), "
                    f"Offset = ({offset_x_g:.4f}, {offset_y_g:.4f}), "
                    f"Magnitude = {total_magnitude_g:.3f}, Physical Offset = {physical_offset_g:.0f} pc"
                )
                print(result_g)
                output_file.write(result_g + "\n")
            else:
                result_g = f"{description}: Failed to extract G (x, y)."
                print(result_g)
                output_file.write(result_g + "\n")

            if x_g2 is not None and y_g2 is not None:
                offset_x_g2, offset_y_g2 = calculate_offset(x_g2, y_g2, observed_g2[0], observed_g2[1])
                total_magnitude_g2 = calculate_total_magnitude(offset_x_g2, offset_y_g2)
                physical_offset_g2 = calculate_physical_offset(total_magnitude_g2, redshift)
                result_g2 = (
                    f"{description} : G2 (x, y) = ({x_g2:.4f}, {y_g2:.4f}), "
                    f"Offset = ({offset_x_g2:.4f}, {offset_y_g2:.4f}), "
                    f"Magnitude = {total_magnitude_g2:.3f}, Physical Offset = {physical_offset_g2:.0f} pc"
                )
                print(result_g2)
                output_file.write(result_g2 + "\n")


# List of file paths to process
file_list = [
    "WFI2033/SP/outSP_optresult.dat",
    "WFI2033/SPG/outSPG_optresult.dat",
    "WFI2033/SPR/outSPR_optresult.dat",
    "WFI2033/SPGR/outSPGR_optresult.dat",
    "WFI2033/SPF/outSPF_optresult.dat",
    "WFI2033/SPFG/outSPFG_optresult.dat",
    "WFI2033/SPFR/outSPFR_optresult.dat",
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
    "WFI2033/NPFGR/outNPFGR_optresult.dat",
]

# Process the files and write results to a text file
process_multiple_files(file_list, "WFI2033_lens_centre_offset_results_gau_priored.txt")


# In[5]:


# SDSSJ1330 lens centre offset code
import math
from astropy.cosmology import Planck18 as cosmo

def extract_coordinates(filename):
    """
    Extract the (x, y) coordinates for the first and second mentioned profiles after the latest "chi^2" line.
    The first profile is labeled as G, and the second profile is labeled as G2.
    """
    x_coordinate_g = None
    y_coordinate_g = None
    x_coordinate_g2 = None
    y_coordinate_g2 = None

    try:
        with open(filename, 'r') as file:
            lines = file.readlines()

        # Find the index of the latest "chi^2" line
        latest_chi2_index = None
        for i, line in enumerate(lines):
            if "chi^2" in line:
                latest_chi2_index = i

        if latest_chi2_index is None:
            return None, None, None, None

        # Process lines after the latest "chi^2" line
        for line in lines[latest_chi2_index + 1:]:
            if "lens   sie" in line or "lens   pow" in line or "lens   anfw" in line:
                # Extract numerical values using a fixed format approach
                parts = [value for value in line.split() if value]

                try:
                    # Extract the x and y coordinates (3rd and 4th parameters)
                    x_coord = float(parts[4])  # 3rd parameter
                    y_coord = float(parts[5])  # 4th parameter

                    # Assign G and G2 based on the order of appearance
                    if x_coordinate_g is None and y_coordinate_g is None:
                        x_coordinate_g = x_coord
                        y_coordinate_g = y_coord
                    elif x_coordinate_g2 is None and y_coordinate_g2 is None:
                        x_coordinate_g2 = x_coord
                        y_coordinate_g2 = y_coord
                        break  # Stop after finding G and G2
                except (ValueError, IndexError):
                    return x_coordinate_g, y_coordinate_g, None, None

    except FileNotFoundError:
        return None, None, None, None
    except Exception:
        return None, None, None, None

    return x_coordinate_g, y_coordinate_g, x_coordinate_g2, y_coordinate_g2


def rename_file(file_path):
    """
    Rename a file based on its naming conventions.
    """
    try:
        parts = file_path.split("/")
        lens_system = parts[0]
        file_name = parts[-1]
        key = file_name.split("_optresult")[0][3:]
        profile_map = {"S": "SIE", "P": "POW", "N": "NFW"}
        profile = profile_map.get(key[0], "Unknown")
        return f"{lens_system} {profile} {key}"
    except IndexError:
        return "Unknown Description"


def calculate_offset(x, y, observed_x, observed_y):
    """
    Calculate the offset in absolute scale.
    """
    return abs(x - observed_x), abs(y - observed_y)


def calculate_total_magnitude(offset_x, offset_y):
    """
    Calculate the total magnitude of the offset using the formula:
    sqrt((x1 - x2)**2 + (y1 - y2)**2).
    """
    return math.sqrt(offset_x**2 + offset_y**2)


def calculate_physical_offset(magnitude_arcsec, redshift):
    """
    Calculate the physical offset in meters using the formula L = r * θ.
    θ is the magnitude in radians, and r is the angular diameter distance.
    """
    theta_radians = math.radians(magnitude_arcsec / 3600.0)
    r_meters = cosmo.angular_diameter_distance(redshift).to('m').value
    l_meters = r_meters * theta_radians
    l_pc = l_meters / 3.086e16
    return l_pc


def calculate_reference_scale(redshift):
    """
    Calculate the arcsecond offset for a physical distance of 10 pc.
    """
    distance_10pc_meters = 10 * 3.086e16
    r_meters = cosmo.angular_diameter_distance(redshift).to('m').value
    theta_radians = distance_10pc_meters / r_meters
    theta_arcseconds = math.degrees(theta_radians) * 3600.0
    return theta_arcseconds


def process_multiple_files(file_list, output_filename):
    """
    Process a list of files, rename them, and write the extracted results including offsets to a file,
    while also printing the results to the console.
    """
    observed_g = (0, 0)
    redshift = 0.373
    reference_scale = calculate_reference_scale(redshift)

    with open(output_filename, 'w') as output_file:
        # Write reference scale at the top
        header = f"Reference Scale: (x, y) = ({reference_scale:.4f}, {reference_scale:.4f}) arcsec for 10 pc\n\n"
        print(header.strip())
        output_file.write(header)

        for file_name in file_list:
            # Rename the file based on naming conventions
            description = rename_file(file_name)

            # Extract x and y coordinates for G and G2
            x_g, y_g, x_g2, y_g2 = extract_coordinates(file_name)

            if x_g is not None and y_g is not None:
                offset_x_g, offset_y_g = calculate_offset(x_g, y_g, observed_g[0], observed_g[1])
                total_magnitude_g = calculate_total_magnitude(offset_x_g, offset_y_g)
                physical_offset_g = calculate_physical_offset(total_magnitude_g, redshift)
                result_g = (
                    f"{description} : G (x, y) = ({x_g:.4f}, {y_g:.4f}), "
                    f"Offset = ({offset_x_g:.4f}, {offset_y_g:.4f}), "
                    f"Magnitude = {total_magnitude_g:.3f}, Physical Offset = {physical_offset_g:.0f} pc"
                )
                print(result_g)
                output_file.write(result_g + "\n")
            else:
                result_g = f"{description}: Failed to extract G (x, y)."
                print(result_g)
                output_file.write(result_g + "\n")

            if x_g2 is not None and y_g2 is not None:
                offset_x_g2, offset_y_g2 = calculate_offset(x_g2, y_g2, observed_g2[0], observed_g2[1])
                total_magnitude_g2 = calculate_total_magnitude(offset_x_g2, offset_y_g2)
                physical_offset_g2 = calculate_physical_offset(total_magnitude_g2, redshift)
                result_g2 = (
                    f"{description} : G2 (x, y) = ({x_g2:.4f}, {y_g2:.4f}), "
                    f"Offset = ({offset_x_g2:.4f}, {offset_y_g2:.4f}), "
                    f"Magnitude = {total_magnitude_g2:.3f}, Physical Offset = {physical_offset_g2:.0f} pc"
                )
                print(result_g2)
                output_file.write(result_g2 + "\n")


# List of file paths to process
file_list = [
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
    "SDSSJ1330/NPFC/outNPFR_optresult.dat"
]

# Process the files and write results to a text file
process_multiple_files(file_list, "SDSSJ1330_lens_centre_offset_results_gau_priored.txt")


# In[6]:


# WFI2026 lens centre offset code
import math
from astropy.cosmology import Planck18 as cosmo

def extract_coordinates(filename):
    """
    Extract the (x, y) coordinates for the first and second mentioned profiles after the latest "chi^2" line.
    The first profile is labeled as G, and the second profile is labeled as G2.
    """
    x_coordinate_g = None
    y_coordinate_g = None
    x_coordinate_g2 = None
    y_coordinate_g2 = None

    try:
        with open(filename, 'r') as file:
            lines = file.readlines()

        # Find the index of the latest "chi^2" line
        latest_chi2_index = None
        for i, line in enumerate(lines):
            if "chi^2" in line:
                latest_chi2_index = i

        if latest_chi2_index is None:
            return None, None, None, None

        # Process lines after the latest "chi^2" line
        for line in lines[latest_chi2_index + 1:]:
            if "lens   sie" in line or "lens   pow" in line or "lens   anfw" in line:
                # Extract numerical values using a fixed format approach
                parts = [value for value in line.split() if value]

                try:
                    # Extract the x and y coordinates (3rd and 4th parameters)
                    x_coord = float(parts[4])  # 3rd parameter
                    y_coord = float(parts[5])  # 4th parameter

                    # Assign G and G2 based on the order of appearance
                    if x_coordinate_g is None and y_coordinate_g is None:
                        x_coordinate_g = x_coord
                        y_coordinate_g = y_coord
                    elif x_coordinate_g2 is None and y_coordinate_g2 is None:
                        x_coordinate_g2 = x_coord
                        y_coordinate_g2 = y_coord
                        break  # Stop after finding G and G2
                except (ValueError, IndexError):
                    return x_coordinate_g, y_coordinate_g, None, None

    except FileNotFoundError:
        return None, None, None, None
    except Exception:
        return None, None, None, None

    return x_coordinate_g, y_coordinate_g, x_coordinate_g2, y_coordinate_g2


def rename_file(file_path):
    """
    Rename a file based on its naming conventions.
    """
    try:
        parts = file_path.split("/")
        lens_system = parts[0]
        file_name = parts[-1]
        key = file_name.split("_optresult")[0][3:]
        profile_map = {"S": "SIE", "P": "POW", "N": "NFW"}
        profile = profile_map.get(key[0], "Unknown")
        return f"{lens_system} {profile} {key}"
    except IndexError:
        return "Unknown Description"


def calculate_offset(x, y, observed_x, observed_y):
    """
    Calculate the offset in absolute scale.
    """
    return abs(x - observed_x), abs(y - observed_y)


def calculate_total_magnitude(offset_x, offset_y):
    """
    Calculate the total magnitude of the offset using the formula:
    sqrt((x1 - x2)**2 + (y1 - y2)**2).
    """
    return math.sqrt(offset_x**2 + offset_y**2)


def calculate_physical_offset(magnitude_arcsec, redshift):
    """
    Calculate the physical offset in meters using the formula L = r * θ.
    θ is the magnitude in radians, and r is the angular diameter distance.
    """
    theta_radians = math.radians(magnitude_arcsec / 3600.0)
    r_meters = cosmo.angular_diameter_distance(redshift).to('m').value
    l_meters = r_meters * theta_radians
    l_pc = l_meters / 3.086e16
    return l_pc


def calculate_reference_scale(redshift):
    """
    Calculate the arcsecond offset for a physical distance of 10 pc.
    """
    distance_10pc_meters = 10 * 3.086e16
    r_meters = cosmo.angular_diameter_distance(redshift).to('m').value
    theta_radians = distance_10pc_meters / r_meters
    theta_arcseconds = math.degrees(theta_radians) * 3600.0
    return theta_arcseconds


def process_multiple_files(file_list, output_filename):
    """
    Process a list of files, rename them, and write the extracted results including offsets to a file,
    while also printing the results to the console.
    """
    observed_g = (0, 0)
    redshift = 1.04
    reference_scale = calculate_reference_scale(redshift)

    with open(output_filename, 'w') as output_file:
        # Write reference scale at the top
        header = f"Reference Scale: (x, y) = ({reference_scale:.4f}, {reference_scale:.4f}) arcsec for 10 pc\n\n"
        print(header.strip())
        output_file.write(header)

        for file_name in file_list:
            # Rename the file based on naming conventions
            description = rename_file(file_name)

            # Extract x and y coordinates for G and G2
            x_g, y_g, x_g2, y_g2 = extract_coordinates(file_name)

            if x_g is not None and y_g is not None:
                offset_x_g, offset_y_g = calculate_offset(x_g, y_g, observed_g[0], observed_g[1])
                total_magnitude_g = calculate_total_magnitude(offset_x_g, offset_y_g)
                physical_offset_g = calculate_physical_offset(total_magnitude_g, redshift)
                result_g = (
                    f"{description} : G (x, y) = ({x_g:.4f}, {y_g:.4f}), "
                    f"Offset = ({offset_x_g:.4f}, {offset_y_g:.4f}), "
                    f"Magnitude = {total_magnitude_g:.3f}, Physical Offset = {physical_offset_g:.0f} pc"
                )
                print(result_g)
                output_file.write(result_g + "\n")
            else:
                result_g = f"{description}: Failed to extract G (x, y)."
                print(result_g)
                output_file.write(result_g + "\n")

            if x_g2 is not None and y_g2 is not None:
                offset_x_g2, offset_y_g2 = calculate_offset(x_g2, y_g2, observed_g2[0], observed_g2[1])
                total_magnitude_g2 = calculate_total_magnitude(offset_x_g2, offset_y_g2)
                physical_offset_g2 = calculate_physical_offset(total_magnitude_g2, redshift)
                result_g2 = (
                    f"{description} : G2 (x, y) = ({x_g2:.4f}, {y_g2:.4f}), "
                    f"Offset = ({offset_x_g2:.4f}, {offset_y_g2:.4f}), "
                    f"Magnitude = {total_magnitude_g2:.3f}, Physical Offset = {physical_offset_g2:.0f} pc"
                )
                print(result_g2)
                output_file.write(result_g2 + "\n")


# List of file paths to process
file_list = [
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
    "WFI2026/NPFC/outNPFR_optresult.dat"
]

# Process the files and write results to a text file
process_multiple_files(file_list, "WFI2026_lens_centre_offset_results_gau_priored.txt")


# In[7]:


# WGDJ0405 lens centre offset code
import math
from astropy.cosmology import Planck18 as cosmo

def extract_coordinates(filename):
    """
    Extract the (x, y) coordinates for the first and second mentioned profiles after the latest "chi^2" line.
    The first profile is labeled as G, and the second profile is labeled as G2.
    """
    x_coordinate_g = None
    y_coordinate_g = None
    x_coordinate_g2 = None
    y_coordinate_g2 = None

    try:
        with open(filename, 'r') as file:
            lines = file.readlines()

        # Find the index of the latest "chi^2" line
        latest_chi2_index = None
        for i, line in enumerate(lines):
            if "chi^2" in line:
                latest_chi2_index = i

        if latest_chi2_index is None:
            return None, None, None, None

        # Process lines after the latest "chi^2" line
        for line in lines[latest_chi2_index + 1:]:
            if "lens   sie" in line or "lens   pow" in line or "lens   anfw" in line:
                # Extract numerical values using a fixed format approach
                parts = [value for value in line.split() if value]

                try:
                    # Extract the x and y coordinates (3rd and 4th parameters)
                    x_coord = float(parts[4])  # 3rd parameter
                    y_coord = float(parts[5])  # 4th parameter

                    # Assign G and G2 based on the order of appearance
                    if x_coordinate_g is None and y_coordinate_g is None:
                        x_coordinate_g = x_coord
                        y_coordinate_g = y_coord
                    elif x_coordinate_g2 is None and y_coordinate_g2 is None:
                        x_coordinate_g2 = x_coord
                        y_coordinate_g2 = y_coord
                        break  # Stop after finding G and G2
                except (ValueError, IndexError):
                    return x_coordinate_g, y_coordinate_g, None, None

    except FileNotFoundError:
        return None, None, None, None
    except Exception:
        return None, None, None, None

    return x_coordinate_g, y_coordinate_g, x_coordinate_g2, y_coordinate_g2


def rename_file(file_path):
    """
    Rename a file based on its naming conventions.
    """
    try:
        parts = file_path.split("/")
        lens_system = parts[0]
        file_name = parts[-1]
        key = file_name.split("_optresult")[0][3:]
        profile_map = {"S": "SIE", "P": "POW", "N": "NFW"}
        profile = profile_map.get(key[0], "Unknown")
        return f"{lens_system} {profile} {key}"
    except IndexError:
        return "Unknown Description"


def calculate_offset(x, y, observed_x, observed_y):
    """
    Calculate the offset in absolute scale.
    """
    return abs(x - observed_x), abs(y - observed_y)


def calculate_total_magnitude(offset_x, offset_y):
    """
    Calculate the total magnitude of the offset using the formula:
    sqrt((x1 - x2)**2 + (y1 - y2)**2).
    """
    return math.sqrt(offset_x**2 + offset_y**2)


def calculate_physical_offset(magnitude_arcsec, redshift):
    """
    Calculate the physical offset in meters using the formula L = r * θ.
    θ is the magnitude in radians, and r is the angular diameter distance.
    """
    theta_radians = math.radians(magnitude_arcsec / 3600.0)
    r_meters = cosmo.angular_diameter_distance(redshift).to('m').value
    l_meters = r_meters * theta_radians
    l_pc = l_meters / 3.086e16
    return l_pc


def calculate_reference_scale(redshift):
    """
    Calculate the arcsecond offset for a physical distance of 10 pc.
    """
    distance_10pc_meters = 10 * 3.086e16
    r_meters = cosmo.angular_diameter_distance(redshift).to('m').value
    theta_radians = distance_10pc_meters / r_meters
    theta_arcseconds = math.degrees(theta_radians) * 3600.0
    return theta_arcseconds


def process_multiple_files(file_list, output_filename):
    """
    Process a list of files, rename them, and write the extracted results including offsets to a file,
    while also printing the results to the console.
    """
    observed_g = (0, 0)
    redshift = 0.29
    reference_scale = calculate_reference_scale(redshift)

    with open(output_filename, 'w') as output_file:
        # Write reference scale at the top
        header = f"Reference Scale: (x, y) = ({reference_scale:.4f}, {reference_scale:.4f}) arcsec for 10 pc\n\n"
        print(header.strip())
        output_file.write(header)

        for file_name in file_list:
            # Rename the file based on naming conventions
            description = rename_file(file_name)

            # Extract x and y coordinates for G and G2
            x_g, y_g, x_g2, y_g2 = extract_coordinates(file_name)

            if x_g is not None and y_g is not None:
                offset_x_g, offset_y_g = calculate_offset(x_g, y_g, observed_g[0], observed_g[1])
                total_magnitude_g = calculate_total_magnitude(offset_x_g, offset_y_g)
                physical_offset_g = calculate_physical_offset(total_magnitude_g, redshift)
                result_g = (
                    f"{description} : G (x, y) = ({x_g:.4f}, {y_g:.4f}), "
                    f"Offset = ({offset_x_g:.4f}, {offset_y_g:.4f}), "
                    f"Magnitude = {total_magnitude_g:.3f}, Physical Offset = {physical_offset_g:.0f} pc"
                )
                print(result_g)
                output_file.write(result_g + "\n")
            else:
                result_g = f"{description}: Failed to extract G (x, y)."
                print(result_g)
                output_file.write(result_g + "\n")

            if x_g2 is not None and y_g2 is not None:
                offset_x_g2, offset_y_g2 = calculate_offset(x_g2, y_g2, observed_g2[0], observed_g2[1])
                total_magnitude_g2 = calculate_total_magnitude(offset_x_g2, offset_y_g2)
                physical_offset_g2 = calculate_physical_offset(total_magnitude_g2, redshift)
                result_g2 = (
                    f"{description} : G2 (x, y) = ({x_g2:.4f}, {y_g2:.4f}), "
                    f"Offset = ({offset_x_g2:.4f}, {offset_y_g2:.4f}), "
                    f"Magnitude = {total_magnitude_g2:.3f}, Physical Offset = {physical_offset_g2:.0f} pc"
                )
                print(result_g2)
                output_file.write(result_g2 + "\n")


# List of file paths to process
file_list = [
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

# Process the files and write results to a text file
process_multiple_files(file_list, "WGDJ0405_lens_centre_offset_results_gau_priored.txt")


# In[9]:


# WGD2038 lens centre offset code
import math
from astropy.cosmology import Planck18 as cosmo

def extract_coordinates(filename):
    """
    Extract the (x, y) coordinates for the first and second mentioned profiles after the latest "chi^2" line.
    The first profile is labeled as G, and the second profile is labeled as G2.
    """
    x_coordinate_g = None
    y_coordinate_g = None
    x_coordinate_g2 = None
    y_coordinate_g2 = None

    try:
        with open(filename, 'r') as file:
            lines = file.readlines()

        # Find the index of the latest "chi^2" line
        latest_chi2_index = None
        for i, line in enumerate(lines):
            if "chi^2" in line:
                latest_chi2_index = i

        if latest_chi2_index is None:
            return None, None, None, None

        # Process lines after the latest "chi^2" line
        for line in lines[latest_chi2_index + 1:]:
            if "lens   sie" in line or "lens   pow" in line or "lens   anfw" in line:
                # Extract numerical values using a fixed format approach
                parts = [value for value in line.split() if value]

                try:
                    # Extract the x and y coordinates (3rd and 4th parameters)
                    x_coord = float(parts[4])  # 3rd parameter
                    y_coord = float(parts[5])  # 4th parameter

                    # Assign G and G2 based on the order of appearance
                    if x_coordinate_g is None and y_coordinate_g is None:
                        x_coordinate_g = x_coord
                        y_coordinate_g = y_coord
                    elif x_coordinate_g2 is None and y_coordinate_g2 is None:
                        x_coordinate_g2 = x_coord
                        y_coordinate_g2 = y_coord
                        break  # Stop after finding G and G2
                except (ValueError, IndexError):
                    return x_coordinate_g, y_coordinate_g, None, None

    except FileNotFoundError:
        return None, None, None, None
    except Exception:
        return None, None, None, None

    return x_coordinate_g, y_coordinate_g, x_coordinate_g2, y_coordinate_g2


def rename_file(file_path):
    """
    Rename a file based on its naming conventions.
    """
    try:
        parts = file_path.split("/")
        lens_system = parts[0]
        file_name = parts[-1]
        key = file_name.split("_optresult")[0][3:]
        profile_map = {"S": "SIE", "P": "POW", "N": "NFW"}
        profile = profile_map.get(key[0], "Unknown")
        return f"{lens_system} {profile} {key}"
    except IndexError:
        return "Unknown Description"


def calculate_offset(x, y, observed_x, observed_y):
    """
    Calculate the offset in absolute scale.
    """
    return abs(x - observed_x), abs(y - observed_y)


def calculate_total_magnitude(offset_x, offset_y):
    """
    Calculate the total magnitude of the offset using the formula:
    sqrt((x1 - x2)**2 + (y1 - y2)**2).
    """
    return math.sqrt(offset_x**2 + offset_y**2)


def calculate_physical_offset(magnitude_arcsec, redshift):
    """
    Calculate the physical offset in meters using the formula L = r * θ.
    θ is the magnitude in radians, and r is the angular diameter distance.
    """
    theta_radians = math.radians(magnitude_arcsec / 3600.0)
    r_meters = cosmo.angular_diameter_distance(redshift).to('m').value
    l_meters = r_meters * theta_radians
    l_pc = l_meters / 3.086e16
    return l_pc


def calculate_reference_scale(redshift):
    """
    Calculate the arcsecond offset for a physical distance of 10 pc.
    """
    distance_10pc_meters = 10 * 3.086e16
    r_meters = cosmo.angular_diameter_distance(redshift).to('m').value
    theta_radians = distance_10pc_meters / r_meters
    theta_arcseconds = math.degrees(theta_radians) * 3600.0
    return theta_arcseconds


def process_multiple_files(file_list, output_filename):
    """
    Process a list of files, rename them, and write the extracted results including offsets to a file,
    while also printing the results to the console.
    """
    observed_g = (0, 0)
    redshift = 0.23
    reference_scale = calculate_reference_scale(redshift)

    with open(output_filename, 'w') as output_file:
        # Write reference scale at the top
        header = f"Reference Scale: (x, y) = ({reference_scale:.4f}, {reference_scale:.4f}) arcsec for 10 pc\n\n"
        print(header.strip())
        output_file.write(header)

        for file_name in file_list:
            # Rename the file based on naming conventions
            description = rename_file(file_name)

            # Extract x and y coordinates for G and G2
            x_g, y_g, x_g2, y_g2 = extract_coordinates(file_name)

            if x_g is not None and y_g is not None:
                offset_x_g, offset_y_g = calculate_offset(x_g, y_g, observed_g[0], observed_g[1])
                total_magnitude_g = calculate_total_magnitude(offset_x_g, offset_y_g)
                physical_offset_g = calculate_physical_offset(total_magnitude_g, redshift)
                result_g = (
                    f"{description} : G (x, y) = ({x_g:.4f}, {y_g:.4f}), "
                    f"Offset = ({offset_x_g:.4f}, {offset_y_g:.4f}), "
                    f"Magnitude = {total_magnitude_g:.3f}, Physical Offset = {physical_offset_g:.0f} pc"
                )
                print(result_g)
                output_file.write(result_g + "\n")
            else:
                result_g = f"{description}: Failed to extract G (x, y)."
                print(result_g)
                output_file.write(result_g + "\n")

            if x_g2 is not None and y_g2 is not None:
                offset_x_g2, offset_y_g2 = calculate_offset(x_g2, y_g2, observed_g2[0], observed_g2[1])
                total_magnitude_g2 = calculate_total_magnitude(offset_x_g2, offset_y_g2)
                physical_offset_g2 = calculate_physical_offset(total_magnitude_g2, redshift)
                result_g2 = (
                    f"{description} : G2 (x, y) = ({x_g2:.4f}, {y_g2:.4f}), "
                    f"Offset = ({offset_x_g2:.4f}, {offset_y_g2:.4f}), "
                    f"Magnitude = {total_magnitude_g2:.3f}, Physical Offset = {physical_offset_g2:.0f} pc"
                )
                print(result_g2)
                output_file.write(result_g2 + "\n")


# List of file paths to process
file_list = [
    "WGD2038/SPC/outSP_optresult.dat",
    "WGD2038/SPC/outSPR_optresult.dat",
    "WGD2038/SPFC/outSPF_optresult.dat",
    "WGD2038/SPFC/outSPFR_optresult.dat",
    "WGD2038/PPC/outPP_optresult.dat",
    "WGD2038/PPC/outPPR_optresult.dat",
    "WGD2038/PPFC/outPPF_optresult.dat",
    "WGD2038/PPFC/outPPFR_optresult.dat",
    "WGD2038/NPC/outNP_optresult.dat",
    "WGD2038/NPC/outNPR_optresult.dat",
    "WGD2038/NPFC/outNPF_optresult.dat",
    "WGD2038/NPFC/outNPFR_optresult.dat"
]

# Process the files and write results to a text file
process_multiple_files(file_list, "WGD2038_lens_centre_offset_results_gau_priored.txt")

