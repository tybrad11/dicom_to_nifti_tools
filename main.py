# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


import convert_all_files_to_nifti
import os

# dicom_path = os.path.join(r"\\onfnas01.uwhis.hosp.wisc.edu\radiology\Research\Bradshaw\Lymphoma_UW_Retrospective\Data\uw_analyzed", 'dicom' )
dicom_path = os.path.join(r"H:\Data\tmp", 'uw_analyzed' )

# nifti_path = os.path.join(r"\\onfnas01.uwhis.hosp.wisc.edu\radiology\Research\Bradshaw\Lymphoma_UW_Retrospective\Data\uw_analyzed", 'nifti')
nifti_path = os.path.join(r"H:\Data\tmp", 'uw_analyzed_nifti')


convert_all_files_to_nifti.convert_all_files_to_nifti(dicom_path, nifti_path, 'PT')


