# -*- coding: utf-8 -*-
"""
Main functions:

- execute_command
- read_spar
- get_rotation_matrice_from_euler_angles
- create_mask_voxel
- convert_to_bids
- placement_new_voxel
"""

import json
import os
import subprocess
from ast import literal_eval

import nibabel as nib
import numpy as np
from scipy.spatial import Delaunay
from scipy.spatial.transform import Rotation


def execute_command(command):
    """
    Execute command

    :param command: command to execute (a list)

    Example:
    - command = ['cd', 'path']
    """
    print("\n", command)
    p = subprocess.Popen(
        command,
        shell=False,
        bufsize=-1,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        close_fds=True,
    )

    print("--------->PID:", p.pid)

    (sdtoutl, stderrl) = p.communicate()
    if str(sdtoutl) != "":
        print("sdtoutl: ", sdtoutl.decode())
    if str(stderrl) != "":
        print("stderrl: ", stderrl.decode())

    result = p.wait()

    return result, stderrl, sdtoutl


def read_spar(filename):
    """Read the .spar file.
    :param filename: file path

    :return: dict of parameters read from spar file
    :rtype: dict
    """

    parameter_dict = {}
    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            # ignore comments (!) and empty lines
            if line == "\n" or line.startswith("!"):
                continue

            # Handle
            key, value = map(str.strip, line.split(":", 1))
            try:
                val = literal_eval(value)
            except (ValueError, SyntaxError):
                if value == "":
                    val = None
                else:
                    val = value

            parameter_dict.update({key: val})

    return parameter_dict


def get_rotation_matrice_from_euler_angles(alpha, beta, gamma, degrees=True):
    """
    Obtains rotation matrice using Euler angles

    :params alpha: first angle (a float)
    :params beta: second angle (a float)
    :params gamma: third angle (a float)
    :params degrees: if not in rad (a boolean)
    :return: rotation matrix
    :rtype: numpy array

    """
    rad = np.pi / 180

    if degrees is True:
        alpha = alpha * rad
        beta = beta * rad
        gamma = gamma * rad

    xrot = np.zeros((3, 3))
    xrot[0, 0] = 1
    xrot[1, 1] = np.cos(alpha)
    xrot[1, 2] = -np.sin(alpha)
    xrot[2, 1] = np.sin(alpha)
    xrot[2, 2] = np.cos(alpha)

    yrot = np.zeros((3, 3))
    yrot[0, 0] = np.cos(beta)
    yrot[0, 2] = np.sin(beta)
    yrot[1, 1] = 1
    yrot[2, 0] = -np.sin(beta)
    yrot[2, 2] = np.cos(beta)

    zrot = np.zeros((3, 3))
    zrot[0, 0] = np.cos(gamma)
    zrot[0, 1] = -np.sin(gamma)
    zrot[1, 0] = np.sin(gamma)
    zrot[1, 1] = np.cos(gamma)
    zrot[2, 2] = 1

    rotation = np.matmul(np.matmul(xrot, yrot), zrot)

    return rotation


def create_mask_voxel(t1_path, parameters, out_path):
    """
    Create a binary mask for a MRS Voxel for Philips
    (using information from .SPAR file)
    Adaptation python of
    https://github.com/markmikkelsen/Gannet/blob/a0753e271069a520c2e6c40c5cee204471ddbb0d/GannetMask_Philips.m#L4

    :param t1_path: file path of the anat data (.nii/ .nii.gz) (a string)
    :param parameters: a dict with parameters from SPAR file. (a dict)
                       Should at least contains info about
                       voxel angulation / center / size
                        ex: parameters = {
                            'ap_size': 16,
                            'lr_size': 16,
                            'cc_size': 16,
                            'lr_off_center': -41.5,
                            'cc_off_center': 9.8,
                            'ap_off_center': -28.27,
                            'ap_angulation': 0.71,
                            'lr_angulation': -1.11,
                            'cc_angulation': 2.08
                        }
    :param out_path: out file path of the binary mask (a string)

    """
    # Read T1
    t1 = nib.load(t1_path)
    dims = t1.shape
    affine = t1.affine
    header = t1.header

    # Generate X, Y, Z coordinate
    X, Y, Z = np.meshgrid(
        np.arange(dims[0]),
        np.arange(dims[1]),
        np.arange(dims[2]),
        indexing="ij",
    )
    XYZ = np.vstack([X.ravel(), Y.ravel(), Z.ravel()])
    # Convert XYZ whith addine matrix
    # XYZ = 3D coordinate of t1
    # XYZ = np.dot(affine[:3, :3], XYZ) + affine[:3, 3:4]
    XYZ = np.matmul(affine[:3, :3], XYZ) + affine[:3, 3:4]

    # In matlab, shifting done, not here
    # % Shift imaging voxel coordinates by half an
    # %imaging voxel so that the XYZ matrix
    # % tells us the x,y,z coordinates of the MIDDLE of that imaging voxel.
    # [~,voxdim2] = spm_get_bbox(V,'fv'); % MM (180220)
    # voxdim2 = abs(voxdim2)';
    # halfpixshift = -voxdim2(1:3)/2;
    # halfpixshift(3) = -halfpixshift(3);
    # XYZ = XYZ + repmat(halfpixshift, [1 size(XYZ,2)]);

    # Get voxel info from spar
    ap_size = parameters["ap_size"]
    lr_size = parameters["lr_size"]
    cc_size = parameters["cc_size"]
    ap_off = parameters["ap_off_center"]
    lr_off = parameters["lr_off_center"]
    cc_off = parameters["cc_off_center"]
    ap_ang = parameters["ap_angulation"]
    lr_ang = parameters["lr_angulation"]
    cc_ang = parameters["cc_angulation"]

    # Flip ap and lr axes to match NIFTI convention
    ap_off = -ap_off
    lr_off = -lr_off
    ap_ang = -ap_ang
    lr_ang = -lr_ang

    vox_ctr = np.array(
        [
            [lr_size / 2, -ap_size / 2, cc_size / 2],
            [-lr_size / 2, -ap_size / 2, cc_size / 2],
            [-lr_size / 2, ap_size / 2, cc_size / 2],
            [lr_size / 2, ap_size / 2, cc_size / 2],
            [-lr_size / 2, ap_size / 2, -cc_size / 2],
            [lr_size / 2, ap_size / 2, -cc_size / 2],
            [lr_size / 2, -ap_size / 2, -cc_size / 2],
            [-lr_size / 2, -ap_size / 2, -cc_size / 2],
        ]
    )

    # Get voxel rotation
    rotation = get_rotation_matrice_from_euler_angles(lr_ang, ap_ang, cc_ang)
    vox_rot = np.matmul(rotation, vox_ctr.T)

    # Get corner coordinate
    vox_ctr_coor = np.array([lr_off, ap_off, cc_off])
    vox_ctr_coor = np.tile(vox_ctr_coor[:, np.newaxis], (1, 8))
    vox_corner = vox_rot + vox_ctr_coor

    # Initial mask and sphere radius
    mask = np.zeros(XYZ.shape[1], dtype=int)
    sphere_radius = np.sqrt(
        (lr_size / 2) ** 2 + (ap_size / 2) ** 2 + (cc_size / 2) ** 2
    )
    dist2voxctr = np.sqrt(
        np.sum(
            (XYZ - np.tile(vox_ctr_coor[:, 0], (XYZ.shape[1], 1)).T) ** 2,
            axis=0,
        )
    )
    sphere_mask = dist2voxctr <= sphere_radius

    mask[sphere_mask] = 1
    XYZ_sphere = XYZ[:, sphere_mask]

    # Triangulation
    tri = Delaunay(np.vstack((vox_corner.T, [lr_off, ap_off, cc_off])))
    tn = tri.find_simplex(XYZ_sphere.T)
    isinside = tn >= 0
    mask[sphere_mask] = isinside

    mask_vol = mask.reshape(dims)
    mask = nib.Nifti2Image(mask_vol, affine, header=header)
    mask.header["descrip"] = "MRS_voxel_mask"
    mask.header["dim"] = t1.header["dim"]
    mask.header["datatype"] = t1.header["datatype"]

    nib.save(mask, out_path)


def convert_to_bids(
    dicom_directory, config_file, out_directory, sub_name, sess_name
):
    """
    Convert to BIDS format (ie convert to NIfTI/json and do the BIDS hierarchy)

    :param dicom_directory: path to dicom (a string)
    :param config_file: path to dcm2bids configuration file (a string)
    :param out_directory: path to out directory (a string)
    :param sub_name: subject name (a string)
    :param sess_name: session name (a string)
    """
    # Launch dcm2bids
    cmd = [
        "dcm2bids",
        "-d",
        dicom_directory,
        "-p",
        sub_name,
        "-s",
        sess_name,
        "-c",
        config_file,
        "-o",
        out_directory,
    ]
    result, stderrl, sdtoutl = execute_command(cmd)


def placement_new_voxel(
    out_directory, t1_ses1_path, list_spar_path, t1_ses2_path
):
    """
    Main function to get the placement of the voxel of the second session.

    This function have been inspired by the following work:
    Woodcock, Eric A. “Automated Voxel Placement:
    A Linux-Based Suite of Tools for Accurate and
    Reliable Single Voxel Coregistration.”
    Journal of Neuroimaging in Psychiatry & Neurology 3, no. 1 (2018): 1–8.
    https://doi.org/10.17756/jnpn.2018-020.

    See code here: https://github.com/ewoodcock/avp_scripts/tree/master


    :param out_directory: path to out directory (a string)
    :param t1_ses1_path: path to T1 NIfTI from session 1 (a string)
    :param list_spar_path: list of the path of the .spar files (a list)
    :param t1_ses2_path: path to T1 NIfTI from session 2 (a string)
    """
    params_new_voxels = []
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

    for i, spar_path in enumerate(list_spar_path):
        out_directory_voxel = os.path.join(
            out_directory, "voxel_" + str(i + 1)
        )
        if not os.path.exists(out_directory_voxel):
            os.makedirs(out_directory_voxel)
        # Create specific out paths by voxel
        mask_path = os.path.join(out_directory_voxel, "voxel_mask.nii")
        mask_ses02_path = os.path.join(
            out_directory_voxel, "voxel_mask_ses02.nii"
        )
        final_info_path = os.path.join(out_directory_voxel, "voxel_ses02.json")

        # Read .SPAR to get session 1 voxel information
        params = read_spar(spar_path)
        create_mask_voxel(t1_ses1_path, params, mask_path)
        # If we want to check the output mask
        # mask_ses02_applyxfm_path = os.path.join(out_directory_voxel,
        # 'voxel_mask_sess02_applyxfm.nii')
        # cmd = ['flirt',  '-in', mask_path,
        # '-ref', t1_ses2_path , '-init',  dof6_inv_matrice, '-out',
        # mask_ses02_applyxfm_path, '-applyxfm']
        # result, stderrl, sdtoutl = execute_command(cmd)

        # Compute final rotation
        subject_rot = np.loadtxt(dof6_inv_matrice)[0:3, 0:3]
        voxel_rotation = get_rotation_matrice_from_euler_angles(
            -params["lr_angulation"],
            -params["ap_angulation"],
            params["cc_angulation"],
        )
        final_rot = np.matmul(voxel_rotation, subject_rot.T)

        # Get new Euler angle
        r = Rotation.from_matrix(final_rot)
        angles = r.as_euler("zyx", degrees=True)
        # It seems that it give the angle as follow : cc, -ap, -lr
        new_cc_ang = angles[0]
        new_ap_ang = angles[1]
        new_lr_ang = angles[2]

        # Get new center coord file
        coord_file = os.path.join(out_directory_voxel, "voxel_coord.txt")
        with open(coord_file, "w", encoding="utf-8") as my_file:
            my_file.write(
                str(-params["lr_off_center"])
                + " "
                + str(-params["ap_off_center"])
                + " "
                + str(params["cc_off_center"])
            )

        cmd = [
            "img2imgcoord",
            "-src",
            t1_ses1_path,
            "-dest",
            t1_ses2_path,
            "-xfm",
            dof6_inv_matrice,
            "-mm",
            coord_file,
        ]
        result, stderrl, sdtoutl = execute_command(cmd)
        new_coord = (
            sdtoutl.decode()
            .replace("Coordinates in Destination volume (in mm)", "")
            .replace("\n", "")
            .split("  ")
        )

        # Create new mask
        params_new = {}
        params_new["ap_size"] = params["ap_size"]
        params_new["lr_size"] = params["lr_size"]
        params_new["cc_size"] = params["cc_size"]
        params_new["lr_off_center"] = round(-float(new_coord[0]), 2)
        params_new["ap_off_center"] = round(-float(new_coord[1]), 2)
        params_new["cc_off_center"] = round(float(new_coord[2]), 2)
        params_new["ap_angulation"] = round(-new_ap_ang, 2)
        params_new["lr_angulation"] = round(-new_lr_ang, 2)
        params_new["cc_angulation"] = round(new_cc_ang, 2)

        create_mask_voxel(t1_ses2_path, params_new, mask_ses02_path)

        # Adapte name for the scanner
        params_new_2 = {}
        params_new_2["AP_size"] = params_new["ap_size"]
        params_new_2["RL_size"] = params_new["lr_size"]
        params_new_2["FH_size"] = params_new["cc_size"]
        params_new_2["RL_off_center"] = params_new["lr_off_center"]
        params_new_2["AP_off_center"] = params_new["ap_off_center"]
        params_new_2["FH_off_center"] = params_new["cc_off_center"]
        params_new_2["AP_angulation"] = params_new["ap_angulation"]
        params_new_2["RL_angulation"] = params_new["lr_angulation"]
        params_new_2["FH_angulation"] = params_new["cc_angulation"]

        with open(final_info_path, "w", encoding="utf-8") as out:
            json.dump(params_new_2, out)

        # Return info to use in the scanner
        params_new_voxels.append(params_new_2)
    return params_new_voxels
