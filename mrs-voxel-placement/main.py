# -*- coding: utf-8 -*-
"""
Module to obtain voxel placement presciption for
a new MRS session using a old MRS session.
"""
import argparse
import json
import os
import sys

from functions_main import placement_new_voxel
from functions_utils import convert_to_bids
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QApplication, QFileDialog, QMainWindow
from PyQt5.uic import loadUi


def main_cli():
    """
    Main CLI
    """
    parser = argparse.ArgumentParser(
        description="Obtain information to position MRS voxel "
        "in longitudinal study"
    )
    parser.add_argument(
        "--session1",
        required=True,
        help="DICOM directory anatomical data session 01",
    )
    parser.add_argument(
        "--session2",
        required=True,
        help="DICOM directory anatomical data session 02",
    )
    parser.add_argument(
        "--spectro_files",
        required=True,
        help="spectro files session 01",
        nargs="+",
    )
    parser.add_argument("--study", required=False, help="Study name")
    parser.add_argument("--patient", required=False, help="Patient name")

    # Set path
    args = parser.parse_args()
    list_spectro_files = args.spectro_files
    t1_ses1_path = args.session1
    t1_ses2_path = args.session2
    study_name = args.study
    patient_name = args.patient

    # Convert DICOM to NIfTI
    config_file = os.path.join(
        os.path.dirname(os.path.realpath(os.path.dirname(__file__))),
        "config",
        "config.json",
    )
    with open(config_file, encoding="utf-8") as my_json:
        data = json.load(my_json)
        bids_config_file = data["BidsConfigFile"]
        out_directory = data["OutputDirectory"]
    out_directory = os.path.join(out_directory, study_name)
    if not os.path.exists(out_directory):
        os.makedirs(out_directory)
    convert_to_bids(
        t1_ses1_path, bids_config_file, out_directory, patient_name, "1"
    )
    convert_to_bids(
        t1_ses2_path, bids_config_file, out_directory, patient_name, "2"
    )

    # Get path NIfTI
    t1_ses1_nifti_path = os.path.join(
        out_directory,
        f"sub-{patient_name}",
        "ses-1",
        "anat",
        f"sub-{patient_name}_ses-1_T1w.nii.gz",
    )
    t1_ses2_nifti_path = os.path.join(
        out_directory,
        f"sub-{patient_name}",
        "ses-2",
        "anat",
        f"sub-{patient_name}_ses-2_T1w.nii.gz",
    )

    # Launch processing to get new voxel placement
    analysis_directory = os.path.join(
        out_directory,
        "derivatives",
        f"sub-{patient_name}",
        "mrs-voxel-placement",
    )
    if not os.path.exists(analysis_directory):
        os.makedirs(analysis_directory)
    manufacturer, params_new_voxels = placement_new_voxel(
        analysis_directory,
        t1_ses1_nifti_path,
        list_spectro_files,
        t1_ses2_nifti_path,
    )

    if "philips" in manufacturer:
        for i, params_new in enumerate(params_new_voxels):
            print(
                f"\n\nFor voxel {str(i + 1)}, the following "
                "information should be used at the scanner: "
                f"\nVoxel Size: \n"
                f' AP: {params_new["AP_size"]:10.2f} \n'
                f' RL: {params_new["RL_size"]:10.2f} \n'
                f' FH: {params_new["FH_size"]:10.2f} '
                f"\nVoxel off center: \n"
                f' AP: {params_new["AP_off_center"]:10.2f} \n'
                f' RL: {params_new["RL_off_center"]:10.2f} \n'
                f' FH: {params_new["FH_off_center"]:10.2f}'
                f"\nVoxel angulation: \n"
                f' AP: {params_new["AP_angulation"]:10.2f} \n'
                f' RL: {params_new["RL_angulation"]:10.2f} \n'
                f' FH: {params_new["FH_angulation"]:10.2f}'
            )
        print("\n\nProcessing done. You can close the application")


class App(QMainWindow):
    """
    Main windows Qt
    """

    def __init__(self):
        super(App, self).__init__()

        # Get ui
        self.dir_code_path = os.path.realpath(os.path.dirname(__file__))
        ui_file = os.path.join(self.dir_code_path, "interface.ui")
        loadUi(ui_file, self)

        config_file = os.path.join(
            os.path.dirname(self.dir_code_path),
            "config",
            "config.json",
        )
        with open(config_file, encoding="utf-8") as my_json:
            data = json.load(my_json)
            self.out_directory = data["OutputDirectory"]
            self.bids_config_file = data["BidsConfigFile"]
        self.spectro_path_1 = ""
        self.spectro_path_2 = ""
        # Connect sigans and slots
        self.pushButton_sess01_T1.clicked.connect(
            lambda: self.browse_directory("sess01_T1", self.out_directory)
        )
        self.pushButton_sess01_spectro_1.clicked.connect(
            lambda: self.browse_file("sess01_spectro_1", self.out_directory)
        )
        self.pushButton_sess01_spectro_2.clicked.connect(
            lambda: self.browse_file("sess01_spectro_2", self.out_directory)
        )
        self.pushButton_sess02_T1.clicked.connect(
            lambda: self.browse_directory("sess02_T1", self.out_directory)
        )
        self.pushButton_run.clicked.connect(self.get_patient_name)
        self.pushButton_run.clicked.connect(self.get_study_name)
        self.pushButton_run.clicked.connect(self.launch_processing)

    def get_patient_name(self):
        """Get patient name from edit line"""
        self.patient_name = self.lineEdit_patient.text()

    def get_study_name(self):
        """Get study name from edit line"""
        self.study_name = self.lineEdit_study.text()

    def browse_file(self, name, out):
        """Browse file"""
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Sélectionner un fichier",
            os.path.dirname(out),
            "Tous les fichiers (*)",
            options=options,
        )
        if file_path:
            if name == "sess01_spectro_1":
                self.spectro_path_1 = file_path
                self.textEdit_sess01_spectro_1.setText(self.spectro_path_1)
            elif name == "sess01_spectro_2":
                self.spectro_path_2 = file_path
                self.textEdit_sess01_spectro_2.setText(self.spectro_path_2)

    def browse_directory(self, name, out):
        """Browse DICOM directory"""
        options = QFileDialog.Options()
        options |= QFileDialog.ShowDirsOnly
        directory = QFileDialog.getExistingDirectory(
            self,
            "Sélectionner un répertoire",
            os.path.dirname(out),
            options=options,
        )

        if directory:
            if name == "sess01_T1":
                self.t1_ses1_path = directory
                self.textEdit_sess01_T1.setText(self.t1_ses1_path)
            elif name == "sess02_T1":
                self.t1_ses2_path = directory
                self.textEdit_sess02_T1.setText(self.t1_ses2_path)

    def launch_processing(self):
        """Launch processing"""
        if self.patient_name and self.study_name:
            if self.spectro_path_1 and self.t1_ses1_path and self.t1_ses2_path:
                # Convert DICOM to NIfTI
                out_directory = os.path.join(
                    self.out_directory, self.study_name
                )
                if not os.path.exists(out_directory):
                    os.makedirs(out_directory)
                convert_to_bids(
                    self.t1_ses1_path,
                    self.bids_config_file,
                    out_directory,
                    self.patient_name,
                    "1",
                )
                convert_to_bids(
                    self.t1_ses2_path,
                    self.bids_config_file,
                    out_directory,
                    self.patient_name,
                    "2",
                )

                # Get path NIfTI
                t1_ses1_nifti_path = os.path.join(
                    out_directory,
                    f"sub-{self.patient_name}",
                    "ses-1",
                    "anat",
                    f"sub-{self.patient_name}_ses-1_T1w.nii.gz",
                )
                t1_ses2_nifti_path = os.path.join(
                    out_directory,
                    f"sub-{self.patient_name}",
                    "ses-2",
                    "anat",
                    f"sub-{self.patient_name}_ses-2_T1w.nii.gz",
                )

                # Launch processing to get new voxel placement
                analysis_directory = os.path.join(
                    out_directory,
                    "derivatives",
                    f"sub-{self.patient_name}",
                    "mrs-voxel-placement",
                )
                list_spectro_files = [self.spectro_path_1]
                if self.spectro_path_2:
                    list_spectro_files.append(self.spectro_path_2)
                if not os.path.exists(analysis_directory):
                    os.makedirs(analysis_directory)
                manufacturer, params_new_voxels = placement_new_voxel(
                    analysis_directory,
                    t1_ses1_nifti_path,
                    list_spectro_files,
                    t1_ses2_nifti_path,
                )

                if "philips" in manufacturer:
                    for i, params_new in enumerate(params_new_voxels):
                        print(
                            f"\n\nFor voxel {str(i + 1)}, the following "
                            "information should be used at the scanner: "
                            f"\nVoxel Size: \n"
                            f' AP: {params_new["AP_size"]:10.2f} \n'
                            f' RL: {params_new["RL_size"]:10.2f} \n'
                            f' FH: {params_new["FH_size"]:10.2f} '
                            f"\nVoxel off center: \n"
                            f' AP: {params_new["AP_off_center"]:10.2f} \n'
                            f' RL: {params_new["RL_off_center"]:10.2f} \n'
                            f' FH: {params_new["FH_off_center"]:10.2f}'
                            f"\nVoxel angulation: \n"
                            f' AP: {params_new["AP_angulation"]:10.2f} \n'
                            f' RL: {params_new["RL_angulation"]:10.2f} \n'
                            f' FH: {params_new["FH_angulation"]:10.2f}'
                        )
                print("\n\nProcessing done. You can close the application")
            else:
                print(
                    "'DICOM T1' fields are mandatory for session 1 and 2. "
                    "'SPAR VOXEL 1' field is mandatory for session 1"
                )
        else:
            print("'Study' and 'Patient identification' fields are mandatory")


def main_gui():
    """Main gui"""
    app = QApplication(sys.argv)
    mainwindow = App()
    widget = QtWidgets.QStackedWidget()
    widget.addWidget(mainwindow)
    widget.setFixedWidth(800)
    widget.setFixedHeight(600)
    widget.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    if len(sys.argv) > 1:
        main_cli()
    else:
        main_gui()
