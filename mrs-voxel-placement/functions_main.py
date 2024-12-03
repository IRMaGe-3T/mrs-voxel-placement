# -*- coding: utf-8 -*-

"""
Main functions:

- coregister_ses02_to_ses01
- placement_new_voxel
"""


import os

from functions_philips import placement_new_voxel_philips
from functions_utils import execute_command, get_value_from_json

# from functions_siemens import placement_new_voxel_siemens


def coregister_ses02_to_ses01(out_directory, t1_ses1_path, t1_ses2_path):
    """
    Registration of ses02 to ses01 using flirt (FSL)
    :param out_directory: path to out directory (a string)
    :param t1_ses1_path: path to T1 NIfTI from session 1 (a string)
    :param t1_ses2_path: path to T1 NIfTI from session 2 (a string)
    """

    # Registration sess02 to sess01 using flirt
    dof6_matrice = os.path.join(out_directory, "ses02_to_ses01_dof6.mat")
    dof6_inv_matrice = dof6_matrice.replace(".mat", "_inverse.mat")
    cmd = [
        "flirt",
        "-dof",
        "6",
        "-in",
        t1_ses2_path,
        "-ref",
        t1_ses1_path,
        "-omat",
        dof6_matrice,
    ]
    result, stderrl, sdtoutl = execute_command(cmd)
    cmd = ["convert_xfm", "-inverse", dof6_matrice, "-omat", dof6_inv_matrice]
    result, stderrl, sdtoutl = execute_command(cmd)

    # In the article, a second affine transformation is computed with dof =9
    # In order to used it with img2imgcoord
    # here it is not done to save time
    # dof9_matrice = os.path.join(os.path.dirname(t1_ses1_path), 'dof9.mat')
    # dof9_inv_matrice = dof9_matrice.replace('.mat', '_inverse.mat')
    # cmd = ['flirt',  '-dof', '9', '-in', t1_ses2_path,
    # '-ref', t1_ses1_path , '-omat',  dof9_matrice]
    # result, stderrl, sdtoutl = execute_command(cmd)
    # cmd = ['convert_xfm', '-inverse', dof9_matrice ,'-omat',
    # dof9_inv_matrice]
    # result, stderrl, sdtoutl = execute_command(cmd)

    return dof6_inv_matrice


def placement_new_voxel(
    out_directory, t1_ses1_path, spectro_files_path, t1_ses2_path
):
    """
    Main function to get the placement of the voxel of the second session.

    This work have been inspired by the following work:
    Woodcock, Eric A. “Automated Voxel Placement:
    A Linux-Based Suite of Tools for Accurate and
    Reliable Single Voxel Coregistration.”
    Journal of Neuroimaging in Psychiatry & Neurology 3, no. 1 (2018): 1–8.
    https://doi.org/10.17756/jnpn.2018-020.
    See code here: https://github.com/ewoodcock/avp_scripts/tree/master

    Modifications of the method habe been done in this code.

    :param out_directory: path to out directory (a string)
    :param t1_ses1_path: path to T1 NIfTI from session 1 (a string)
    :param spectro_files_path: list of the path of the spectro files (a list)
    :param t1_ses2_path: path to T1 NIfTI from session 2 (a string)
    :return:  manufacturer, params_new_voxels
    """
    params_new_voxels = []

    # Check Manufacturer (only Siemens or Philips ok for now)
    manufacturer = get_value_from_json(
        t1_ses1_path.replace(".nii.gz", ".json"), "Manufacturer"
    ).lower()

    # if "philips" not in manufacturer and "siemens" not in manufacturer:
    if "philips" not in manufacturer:
        print(f"Manufactuer {manufacturer} not implemented for this code")
        return
    for i, spectro_path in enumerate(spectro_files_path):
        # Check spectro files format
        if "philips" in manufacturer:
            # Check spar file
            ext = spectro_path.split(".")[-1]
            if ext not in ["SPAR", "spar"]:
                print(
                    "For Philips system, "
                    "spectro files should be '.SPAR' files"
                )
                return
        # elif "siemens" in manufacturer:
        #     ds = pydicom.read_file(spectro_path)
        #     # Check SOP class Storage
        #     storage_method = ds[Tag(0x08, 0x16)].value[0]
        #     # Check number
        #     if storage_method not in
        #       ["MR Spectroscopy Storage", "1.2.840.10008.5.1.4.1.1.4.2"]:
        #         print(
        #             "For Siemens system, "
        #             "spectro files should be DICOM MR Spectroscopy Storage"
        #         )
        #         return

    # Coregister ses02 to ses01 (using flirt FSL)
    dof6_inv_matrice = coregister_ses02_to_ses01(
        out_directory, t1_ses1_path, t1_ses2_path
    )

    for i, spectro_path in enumerate(spectro_files_path):
        # For each spectro file get info for session 2
        out_directory_voxel = os.path.join(
            out_directory, "voxel_" + str(i + 1)
        )
        if not os.path.exists(out_directory_voxel):
            os.makedirs(out_directory_voxel)

        params_new_voxel = {}

        if "philips" in manufacturer:
            # Get info to use in the scanner
            params_new_voxel = placement_new_voxel_philips(
                spectro_path,
                t1_ses1_path,
                t1_ses2_path,
                dof6_inv_matrice,
                out_directory_voxel,
            )
        # elif "siemens" in manufacturer:
        #     # Get info to use in the scanner
        #     params_new_voxel = placement_new_voxel_siemens(
        #         spectro_path,
        #         t1_ses1_path,
        #         t1_ses2_path,
        #         dof6_inv_matrice,
        #         out_directory_voxel
        #     )

        params_new_voxels.append(params_new_voxel)

    return manufacturer, params_new_voxels
