# -*- coding: utf-8 -*-

import json
import math
import os

import nibabel as nib
import numpy as np
import pydicom
from functions_utils import execute_command
from pydicom.tag import Tag
from scipy.spatial import Delaunay


def create_mask_voxel_siemens(t1_path, parameters, out_path):
    """
    Create a binary mask for a MRS Voxel for Siemens
    (using information from DICOM file  MRSpectroscopyStrorage)
    Adaptation python of
    https://github.com/markmikkelsen/Gannet/blob/a0753e271069a520c2e6c40c5cee204471ddbb0d/GannetMask_Philips.m#L4

    :param t1_path: file path of the anat data (.nii/ .nii.gz) (a string)
    :param parameters: a dict with parameters from SPAR file. (a dict)
                       Should at least contains info about
                       voxel angulation / center / size
                        ex: parameters = {
                            "ap_size": 16,
                            "lr_size": 16,
                            "size3": 16,
                            "lr_off_center": -41.5,
                            "cc_off_center": 9.8,
                            "ap_off_center": -28.27,
                            "ap_angulation": 0.71,
                            "lr_angulation": -1.11,
                            "cc_angulation": 2.08
                        }
    :param out_path: out file path of the binary mask (a string)

    """
    # /!\ WIP
    # Read T1
    t1 = nib.load(t1_path)
    dims = t1.shape
    affine = t1.affine
    header = t1.header

    # Get voxel info
    size_thickness = parameters["size_thickness"]
    size_pe_fov = parameters["size_pe_fov"]
    size_ro_fov = parameters["size_ro_fov"]
    # TODO: check order for vox dim
    voxel_dim = [size_pe_fov, size_ro_fov, size_thickness]
    # image_orientation = parameters["image_orientation"]
    # image_position = parameters["image_position"]
    center_lr = parameters["position_center_sag"]
    center_pa = parameters["position_center_cor"]
    center_fh = parameters["position_center_tra"]
    norm_sag = parameters["norm_sag"]
    norm_cor = parameters["norm_cor"]
    norm_tra = parameters["norm_tra"]
    rotation = parameters["rotation"]

    # Get vox orient and get rotation matrix
    # TODO: - or not - ??
    norm = [-norm_sag, -norm_cor, norm_tra]
    # norm = [norm_sag, norm_cor, norm_tra]
    rotation, rotation_matrix = get_rotation_matrix_siemens(norm, rotation)

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
    # [~,voxdim2] = spm_get_bbox(V,"fv"); % MM (180220)
    # voxdim2 = abs(voxdim2)";
    # halfpixshift = -voxdim2(1:3)/2;
    # halfpixshift(3) = -halfpixshift(3);
    # XYZ = XYZ + repmat(halfpixshift, [1 size(XYZ,2)]);

    # Center position given in LPH+
    # Flip ap and lr axes to match NIFTI convention
    center_lr = -center_lr
    center_pa = -center_pa

    vox_ctr = np.array(
        [
            [voxel_dim[0] / 2, -voxel_dim[1] / 2, voxel_dim[2] / 2],
            [-voxel_dim[0] / 2, -voxel_dim[1] / 2, voxel_dim[2] / 2],
            [-voxel_dim[0] / 2, voxel_dim[1] / 2, voxel_dim[2] / 2],
            [voxel_dim[0] / 2, voxel_dim[1] / 2, voxel_dim[2] / 2],
            [-voxel_dim[0] / 2, voxel_dim[1] / 2, -voxel_dim[2] / 2],
            [voxel_dim[0] / 2, voxel_dim[1] / 2, -voxel_dim[2] / 2],
            [voxel_dim[0] / 2, -voxel_dim[1] / 2, -voxel_dim[2] / 2],
            [-voxel_dim[0] / 2, -voxel_dim[1] / 2, -voxel_dim[2] / 2],
        ]
    )

    # Get voxel rotation
    vox_rot = np.matmul(rotation_matrix, vox_ctr.T)

    # Get corner coordinate
    vox_ctr_coor = np.array([center_lr, center_pa, center_fh])
    vox_ctr_coor = np.tile(vox_ctr_coor[:, np.newaxis], (1, 8))
    vox_corner = vox_rot + vox_ctr_coor

    # Initial mask and sphere radius
    mask = np.zeros(XYZ.shape[1], dtype=int)
    sphere_radius = np.sqrt(
        (voxel_dim[0] / 2) ** 2
        + (voxel_dim[1] / 2) ** 2
        + (voxel_dim[2] / 2) ** 2
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
    tri = Delaunay(
        np.vstack((vox_corner.T, [center_lr, center_pa, center_fh]))
    )
    tn = tri.find_simplex(XYZ_sphere.T)
    isinside = tn >= 0
    mask[sphere_mask] = isinside

    mask_vol = mask.reshape(dims)
    mask = nib.Nifti2Image(mask_vol, affine, header=header)
    mask.header["descrip"] = "MRS_voxel_mask"
    mask.header["dim"] = t1.header["dim"]
    mask.header["datatype"] = t1.header["datatype"]

    nib.save(mask, out_path)


def get_rotation_matrix_siemens(norm, rotation=None, norm_2=None):
    # /!\ WIP
    # Find largest element of normal vector to determine primary orientaion
    # Compute phase reference vector
    norm_abs = [abs(norm[0]), abs(norm[1]), abs(norm[2])]
    index_max_norm = norm_abs.index(max(norm_abs))
    phase = [0, 0, 0]

    if index_max_norm == 2:
        # For transversal voxel orientation,
        # the phase reference vector lies in the sagittal plane
        vox_orient = "t"
        phase = [
            0,
            norm[2] * (1 / (norm[1] ** 2 + norm[2] ** 2)) ** 0.5,
            -norm[1] * (1 / (norm[1] ** 2 + norm[2] ** 2)) ** 0.5,
        ]
    elif index_max_norm == 1:
        # For coronal voxel orientation,
        # the phase reference vector lies in the transversal plane
        vox_orient = "c"
        phase = [
            norm[1] * (1 / (norm[0] ** 2 + norm[1] ** 2)) ** 0.5,
            -norm[0] * (1 / (norm[0] ** 2 + norm[1] ** 2)) ** 0.5,
            0,
        ]
    elif index_max_norm == 0:
        # For sagittal voxel orientation,
        # the phase reference vector lies in the transversal plane
        vox_orient = "s"
        phase = [
            -norm[1] * (1 / (norm[0] ** 2 + norm[1] ** 2)) ** 0.5,
            norm[0] * (1 / (norm[0] ** 2 + norm[1] ** 2)) ** 0.5,
            0,
        ]
    print("Vox orient", vox_orient)
    if not rotation and norm_2:
        # Compute rotation
        cross_product = np.cross(phase, norm_2)
        dot_product = np.dot(cross_product, norm)
        if dot_product <= 0:
            rotation = np.arccos(np.dot(norm_2, phase))
            # rotation = -np.arccos(np.dot(norm_2, phase))
        else:
            rotation = -np.arccos(np.dot(norm_2, phase))
            # rotation = np.arccos(np.dot(norm_2, phase))

    # The readout reference vector is the cross product of Norm and Phase
    readout = np.cross(norm, phase)
    # Compute total rotation matrix
    mat1 = np.zeros((4, 4))
    mat1[0:3, 0] = phase
    mat1[0:3, 1] = readout
    mat1[0:3, 2] = norm

    mat_rot = np.array(
        [
            [np.cos(rotation), np.sin(rotation), 0],
            [-np.sin(rotation), np.cos(rotation), 0],
            [0, 0, 1],
        ]
    )

    mat2 = np.dot(mat1[0:3, 0:3], mat_rot)
    mat1[0:3, 0:3] = mat2
    # The MGH vox2ras matrix inverts the Readout column
    mat1 = np.dot(
        mat1,
        np.array([[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
    )

    rotation_matrix = mat1[0:3, 0:3]

    return rotation, rotation_matrix


def read_siemens_dicom(dicom_path):
    """
    Read DICOM spectro Siemens
    """

    params = {}
    ds = pydicom.read_file(dicom_path)

    # For Siemens VE software:
    # the orientation of the voxel is defined by a by a normal vector
    # of the voxel and a rotation around it
    # For Siemens XA software:
    # the orientation of the voxel is defined by a three slabs
    # of excitation/refocusing, so we need to compute orientation

    storage_method = ds[Tag(0x08, 0x16)].value

    if storage_method in [
        "MR Spectroscopy Storage",
        "1.2.840.10008.5.1.4.1.1.4.2",
    ]:
        params["size_thickness"] = ds[Tag(0x018, 0x9126)][0][
            Tag(0x018, 0x9104)
        ].value
        params["size_pe_fov"] = ds[Tag(0x018, 0x9126)][1][
            Tag(0x018, 0x9104)
        ].value
        params["size_ro_fov"] = ds[Tag(0x018, 0x9126)][2][
            Tag(0x018, 0x9104)
        ].value

        slab_orientation_1 = ds[Tag(0x018, 0x9126)][0][
            Tag(0x018, 0x9105)
        ].value
        slab_orientation_2 = ds[Tag(0x018, 0x9126)][1][
            Tag(0x018, 0x9105)
        ].value
        slab_orientation_3 = ds[Tag(0x018, 0x9126)][2][
            Tag(0x018, 0x9105)
        ].value

        params["norm_sag"] = slab_orientation_1[0]
        params["norm_cor"] = slab_orientation_1[1]
        params["norm_tra"] = slab_orientation_1[2]

        params["slab_orientation_1"] = slab_orientation_1
        params["slab_orientation_2"] = slab_orientation_2
        params["slab_orientation_3"] = slab_orientation_3

        mid_slab_position = ds[Tag(0x018, 0x9126)][0][Tag(0x018, 0x9106)].value
        params["position_center_sag"] = mid_slab_position[0]
        params["position_center_cor"] = mid_slab_position[1]
        params["position_center_tra"] = mid_slab_position[2]

        params["image_orientation"] = ds[Tag(0x5200, 0x9230)][0][
            Tag(0x020, 0x9116)
        ][0][Tag(0x020, 0x0037)].value
        params["image_position"] = ds[Tag(0x5200, 0x9230)][0][
            Tag(0x020, 0x9113)
        ][0][Tag(0x020, 0x0032)].value

        # Get roration
        norm = slab_orientation_1
        # norm = [-slab_orientation_1[0],
        # -slab_orientation_1[1], slab_orientation_1[2]]
        norm2 = slab_orientation_2
        # norm2 = [-slab_orientation_2[0],
        # -slab_orientation_2[1], slab_orientation_2[2]]
        rotation, rotation_matrix = get_rotation_matrix_siemens(
            norm, norm_2=norm2
        )
        params["rotation"] = rotation
        rotation_deg = rotation * 180 / math.pi
        params["rotation_deg"] = rotation_deg

    # elif storage_method in
    # TODO: add VE system

    return params


def placement_new_voxel_siemens(
    spectro_path,
    t1_ses1_path,
    t1_ses2_path,
    dof6_inv_matrice,
    out_directory_voxel,
):
    """
    Get information about localization of the voxel for session 2
    for Siemens system
    """
    # /!\ WIP
    # Create specific out paths by voxel
    mask_ses01_path = os.path.join(out_directory_voxel, "voxel_mask_ses01.nii")
    mask_ses02_path = os.path.join(out_directory_voxel, "voxel_mask_ses02.nii")
    final_info_path = os.path.join(
        out_directory_voxel, "info_voxel_final_ses02.json"
    )

    # Read DICOM to get session 1 voxel information
    params = read_siemens_dicom(spectro_path)
    create_mask_voxel_siemens(t1_ses1_path, params, mask_ses01_path)
    # If we want to check the output mask
    # mask_ses02_applyxfm_path = os.path.join(out_directory_voxel,
    # "voxel_mask_sess02_applyxfm.nii")
    # cmd = ["flirt",  "-in", mask_ses01_path,
    # "-ref", t1_ses2_path , "-init",  dof6_inv_matrice, "-out",
    # mask_ses02_applyxfm_path, "-applyxfm"]
    # result, stderrl, sdtoutl = execute_command(cmd)

    # Compute final rotation
    subject_rot = np.loadtxt(dof6_inv_matrice)[0:3, 0:3]

    norm_sag = params["norm_sag"]
    norm_cor = params["norm_cor"]
    norm_tra = params["norm_tra"]
    # TODO: - or not - ??
    norm = [-norm_sag, -norm_cor, norm_tra]
    # norm = [norm_sag, norm_cor, norm_tra]
    rotation = params["rotation"]
    rotation, voxel_rotation = get_rotation_matrix_siemens(norm, rotation)
    final_rot = np.matmul(voxel_rotation, subject_rot.T)
    print("Final rotation", final_rot)

    # Get new rotation and norm vect ??

    # Get new center coord file
    coord_file = os.path.join(out_directory_voxel, "voxel_coord.txt")
    with open(coord_file, "w", encoding="utf-8") as my_file:
        my_file.write(
            str(params["position_center_sag"])
            + " "
            + str(params["position_center_cor"])
            + " "
            + str(params["position_center_tra"])
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
    os.remove(coord_file)

    # # Create new mask
    params_new = {}
    params_new["size_thickness"] = params["size_thickness"]
    params_new["size_pe_fov"] = params["size_pe_fov"]
    params_new["size_ro_fov"] = params["size_ro_fov"]
    params_new["position_center_sag"] = new_coord[0]
    params_new["position_center_cor"] = new_coord[1]
    params_new["position_center_tra"] = new_coord[2]
    # params_new["norm_sag"] =
    # params_new["norm_cor"] =
    # params_new["norm_tra"] =
    # params_new["rotation"] =
    # params_new["rotation_deg"]=

    create_mask_voxel_siemens(t1_ses1_path, params_new, mask_ses02_path)

    # Adapte name for the Siemens scanner
    params_new_2 = {}

    with open(final_info_path, "w", encoding="utf-8") as out:
        json.dump(params_new_2, out)

    return params_new_2
