'''
Module to obtain voxel placement presciption for 
a new MRS session using a old MRS session.
'''
import argparse
import json
import os
import sys

from PyQt5 import QtWidgets
from PyQt5.QtCore import QDir
from PyQt5.QtWidgets import QApplication, QFileDialog, QMainWindow, QMessageBox
from PyQt5.uic import loadUi

from functions import convert_to_bids, placement_new_voxel


def main_cli():
    ''' 
    Main CLI
    '''
    parser = argparse.ArgumentParser(
        description='Obtain information to position MRS voxel '
        'in longitudinal study'
    )
    parser.add_argument(
        '--session1', required=True, help='DICOM directory anatomical data session 01'
    )
    parser.add_argument(
        '--session2', required=True, help='DICOM directory anatomical data session 02'
    )
    parser.add_argument(
        '--spar', required=True, help='spar file session 01', nargs='+',

    )
    parser.add_argument(
        '--study', required=False, help='Study name'
    )
    parser.add_argument(
        '--patient', required=False, help='Patient name'
    )

    # Set path
    args = parser.parse_args()
    list_spar_path = args.spar
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
        t1_ses1_path, bids_config_file, out_directory, patient_name, "1")
    convert_to_bids(
        t1_ses2_path, bids_config_file, out_directory, patient_name, "2")

    # Get path NIfTI
    t1_ses1_nifti_path = os.path.join(
        out_directory,
        f'sub-{patient_name}',
        'ses-1',
        'anat',
        f'sub-{patient_name}_ses-1_T1w.nii.gz'
    )
    t1_ses2_nifti_path = os.path.join(
        out_directory,
        f'sub-{patient_name}',
        'ses-2',
        'anat',
        f'sub-{patient_name}_ses-2_T1w.nii.gz'
    )


    # Launch processing to get new voxel placement
    analysis_directory = os.path.join(
        out_directory,
        'derivatives',
        f'sub-{patient_name}',
        'mrs-voxel-placement'     

    )
    if not os.path.exists(analysis_directory):
        os.makedirs(analysis_directory)
    placement_new_voxel(analysis_directory, t1_ses1_nifti_path,
                        list_spar_path, t1_ses2_nifti_path)


class App(QMainWindow):
    '''
    Main windows Qt
    '''

    def __init__(self):
        super(App, self).__init__()

        # Get ui
        self.dir_code_path = os.path.realpath(os.path.dirname(__file__))
        ui_file = os.path.join(self.dir_code_path, 'interface.ui')
        loadUi(ui_file, self)

        # Connect sigans and slots
        self.pushButton_sess01_T1.clicked.connect(
            lambda: self.browse_directory('sess01_T1'))
        self.pushButton_sess01_spar_1.clicked.connect(
            lambda: self.browse_file('sess01_spar_1'))
        self.pushButton_sess01_spar_2.clicked.connect(
            lambda: self.browse_file('sess01_spar_2'))
        self.pushButton_sess02_T1.clicked.connect(
            lambda: self.browse_directory('sess02_T1'))
        self.pushButton_run.clicked.connect(self.get_patient_name)
        self.pushButton_run.clicked.connect(self.get_study_name)
        self.pushButton_run.clicked.connect(self.launch_processing)

    def get_patient_name(self):
        '''Get patient name from edit line '''
        self.patient_name = self.lineEdit_patient.text()

    def get_study_name(self):
        '''Get study name from edit line '''
        self.study_name = self.lineEdit_study.text()

    def browse_file(self, name):
        '''Browse file'''
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly
        file_path, _ = QFileDialog.getOpenFileName(
            self, 'Sélectionner un fichier',
            QDir.homePath(), 'Tous les fichiers (*)', options=options)

        if file_path:
            if name == 'sess01_spar_1':
                self.spar_path_1 = file_path
                self.textEdit_sess01_spar_1.setText(self.spar_path_1)
            elif name == 'sess01_spar_2':
                self.spar_path_2 = file_path
                self.textEdit_sess01_spar_2.setText(self.spar_path_2)

    def browse_directory(self, name):
        '''Browse DICOM directory'''
        directory = QFileDialog.getExistingDirectory(
            self, 'Sélectionner un répertoire', QDir.homePath())

        if directory:
            if name == 'sess01_T1':
                self.t1_ses1_path = directory
                self.textEdit_sess01_T1.setText(self.t1_ses1_path)
            elif name == 'sess02_T1':
                self.t1_ses2_path = directory
                self.textEdit_sess02_T1.setText(self.t1_ses2_path)

    def launch_processing(self):
        '''Launch processing'''
        if self.patient_name and self.study_name:
            if self.spar_path_1 and self.t1_ses1_path and self.t1_ses2_path:

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
                out_directory = os.path.join(out_directory, self.study_name)
                if not os.path.exists(out_directory):
                    os.makedirs(out_directory)
                convert_to_bids(
                    self.t1_ses1_path, bids_config_file, out_directory,  self.patient_name, "1")
                convert_to_bids(
                    self.t1_ses2_path, bids_config_file, out_directory,  self.patient_name, "2")

                # Get path NIfTI
                t1_ses1_nifti_path = os.path.join(
                    out_directory,
                    f'sub-{self.patient_name}',
                    'ses-1',
                    'anat',
                    f'sub-{self.patient_name}_ses-1_T1w.nii.gz'
                )
                t1_ses2_nifti_path = os.path.join(
                    out_directory,
                    f'sub-{self.patient_name}',
                    'ses-2',
                    'anat',
                    f'sub-{self.patient_name}_ses-2_T1w.nii.gz'
                )

                # Launch processing to get new voxel placement
                analysis_directory = os.path.join(
                    out_directory,
                    'derivatives',
                    f'sub-{self.patient_name}',
                    'mrs-voxel-placement'

                )
                list_spar_files = [self.spar_path_1]
                if self.spar_path_2:
                    list_spar_files.append(self.spar_path_2)
                if not os.path.exists(analysis_directory):
                    os.makedirs(analysis_directory)
                params_new_voxels = placement_new_voxel(analysis_directory, t1_ses1_nifti_path,
                                                        list_spar_files, t1_ses2_nifti_path)

                for i, params_new in enumerate(params_new_voxels):
                    print(
                        f'\n\nFor voxel {str(i + 1)}, '
                        'the following information should be used at the scanner: '
                        f'\nVoxel Size: \n  ap_size : {params_new["ap_size"]} \n  '
                        f'lr_size: {params_new["lr_size"]} \n  cc_size: {params_new["cc_size"]} '
                        f'\nVoxel off center: \n  ap_off_center : {params_new["ap_off_center"]} '
                        f'\n  lr_off_center: {params_new["lr_off_center"]} \n  '
                        f'cc_off_center: {params_new["cc_off_center"]}'
                        f'\nVoxel angulation: \n  ap_angulation : {params_new["ap_angulation"]} '
                        f'\n  lr_angulation: {params_new["lr_angulation"]} \n  '
                        f'cc_angulation: {params_new["cc_angulation"]}'
                    )

def main_gui():
    '''Main gui '''
    app = QApplication(sys.argv)
    mainwindow = App()
    widget = QtWidgets.QStackedWidget()
    widget.addWidget(mainwindow)
    widget.setFixedWidth(800)
    widget.setFixedHeight(600)
    widget.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    if len(sys.argv) > 1:
        main_cli()
    else:
        main_gui()
