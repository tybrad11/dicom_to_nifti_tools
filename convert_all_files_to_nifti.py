
import os
import pydicom
import dicom2nifti
from platipy.dicom.io import rtstruct_to_nifti
import get_all_terminal_subfolders

def isdicom(file_path):
    #borrowed from pydicom filereader.py
    with open(file_path, 'rb') as fp:
        preamble = fp.read(128)
        magic = fp.read(4)
        if magic != b"DICM":
            return False
        else:
            return True

def find_path_to_dicom_image_that_corresponds_with_rtsrtuct(dir_i, dicom_study_date, modality_of_interest, subdirs):
    one_folder_up = os.path.dirname(dir_i)
    patient_folders = [s for s in subdirs if one_folder_up in s]
    corresponding_dicom = ''
    for dir_p in patient_folders:
        files2 = os.listdir(dir_p)
        temp_file = files2[0]
        if isdicom(os.path.join(dir_p, temp_file)) == False:
            continue
        test2_dicom = pydicom.dcmread(os.path.join(dir_p, temp_file))
        if  test2_dicom['00080020'].value == dicom_study_date and test2_dicom['00080060'].value == modality_of_interest:
            corresponding_dicom = dir_p
    return corresponding_dicom



def convert_all_files_to_nifti(top_dicom_folder, top_nifti_folder, modality_of_interest):
    subdirs = get_all_terminal_subfolders(top_dicom_folder)

    for dir_i in subdirs:
        files = os.listdir(dir_i)
        #REMOVE MIM CACHE FILESc
        if len(files) < 1:
            print('Empty folder: '+dir_i)
            continue
        if isdicom(os.path.join(dir_i, files[0])) == False:
            print('Non-DICOM: ' +dir_i)

        #get dicom info for saving
        test_dicom = pydicom.dcmread(os.path.join(dir_i, files[0]))
        dicom_modality = test_dicom['00080060'].value
        dicom_name = str(test_dicom['00100010'].value).lower()
        dicom_id = test_dicom['00100020'].value.lower()
        dicom_series_date = test_dicom['00080021'].value
        dicom_series_description = test_dicom['0008103e'].value

        #unique names for subjects and scans
        subject_save_name = dicom_id + '_' + dicom_name.replace(' ', '_').replace('__', '_')
        subject_save_folder = os.path.join(top_nifti_folder, subject_save_name)
        scan_save_name = dicom_series_date + '_' + dicom_modality + '_' + dicom_series_description.replace(' ', '_')


        if not os.path.exists(subject_save_folder):
            os.makedirs(subject_save_folder)

        if dicom_modality in ['PT', 'CT', 'MR', 'NM']:
            dicom2nifti.dicom_series_to_nifti(dir_i, os.path.join(subject_save_folder, scan_save_name + '.nii.gz'), reorient_nifti=False)

        elif dicom_modality == 'RTSTRUCT':
            #might be multiple rtstructs in folder
            for file_i in files:
                if isdicom(os.path.join(dir_i, file_i)) == True:
                    rt_dicom = pydicom.dcmread(os.path.join(dir_i, file_i))
                    if rt_dicom['00080060'].value == 'RTSTRUCT':
                        dicom_modality = rt_dicom['00080060'].value
                        dicom_name = str(rt_dicom['00100010'].value).lower()
                        dicom_id = rt_dicom['00100020'].value.lower()
                        dicom_study_date = rt_dicom['00080020'].value
                        dicom_series_description = rt_dicom['0008103e'].value

                        # unique names for subjects and scans
                        subject_save_name = dicom_id + '_' + dicom_name.replace(' ', '_').replace('__', '_')
                        subject_save_folder = os.path.join(top_nifti_folder, subject_save_name)
                        scan_save_name = dicom_study_date + '_' + dicom_modality + '_' + dicom_series_description.replace(' ','_')

                        #find the corresponding DICOM series
                        corresp_dicom_path = find_path_to_dicom_image_that_corresponds_with_rtsrtuct(dir_i , dicom_study_date,  modality_of_interest, subdirs)
                        #save
                        rtstruct_nifti_save_path = os.path.join(subject_save_folder, scan_save_name)
                        if not os.path.exists(rtstruct_nifti_save_path):
                            os.makedirs(rtstruct_nifti_save_path)
                        rtstruct_to_nifti.convert_rtstruct(corresp_dicom_path, os.path.join(dir_i, file_i), output_dir = rtstruct_nifti_save_path)
        else:
            print('Modality is not standard modality!')

