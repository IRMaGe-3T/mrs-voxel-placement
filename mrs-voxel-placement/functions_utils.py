# -*- coding: utf-8 -*-
"""
Functions:

- convert_to_bids
- execute_command
- get_euler_angles_fom_rotation_matrix
- get_rotation_matrice_from_euler_angles
- get_value_from_json

"""

import json
import subprocess

import numpy as np


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


def execute_command(command):
    """
    Execute command

    :param command: command to execute (a list)
    :return: result, stderrl, sdtoutl
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


def get_euler_angles_fom_rotation_matrix(mat_rot):
    """
    Convert a rotation matrix in Euler angle(yaw, pitch, roll) (degrees)

    :params mat_rot : rotation matrix 3x3
    :return: a tuple (yaw, pitch, roll) (degrees).
    """
    # Get pitch yaw roll
    if np.isclose(mat_rot[2, 0], -1.0):
        pitch = np.pi / 2
        yaw = np.arctan2(mat_rot[0, 1], mat_rot[0, 2])
        roll = 0
    elif np.isclose(mat_rot[2, 0], 1.0):
        pitch = -np.pi / 2
        yaw = np.arctan2(-mat_rot[0, 1], -mat_rot[0, 2])
        roll = 0
    else:
        # Compute Euler Angle with convention ZYX
        pitch = np.arcsin(-mat_rot[2, 0])
        yaw = np.arctan2(mat_rot[2, 1], mat_rot[2, 2])
        roll = np.arctan2(mat_rot[1, 0], mat_rot[0, 0])

    # Radian to degrees
    yaw = np.degrees(yaw)
    pitch = np.degrees(pitch)
    roll = np.degrees(roll)

    return [yaw, pitch, roll]


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


def get_value_from_json(json_path, tag):
    """
    Get value from a json

    :param json_path: path to json (a string)
    :param tag (a string)
    :return: value
    """

    with open(json_path, "r", encoding="utf-8") as json_file:
        data = json.load(json_file)

        if tag in list(data.keys()):
            value = data[tag]
        else:
            value = "NotFound"
    return value
