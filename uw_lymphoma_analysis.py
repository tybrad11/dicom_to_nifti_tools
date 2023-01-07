from convert_all_files_to_nifti import convert_all_files_to_nifti, convert_specific_rtstruct_to_nifti
import os
from get_all_terminal_subfolders import get_all_terminal_subfolders
from nifti_roi_tools import combine_rois_into_one, combine_rois_into_groups, resample_ct_nii_to_match_pet_nii
import nibabel as nib
import numpy as np
from list_converted_nifti_files import list_converted_nifti_files, find_which_rtstructs_werent_converted
import glob




def one_time_add_patientid_substring_to_nifti_filenames(parent_folder):
    contents = os.listdir(parent_folder)

    for patient in contents:
        if os.path.isdir(os.path.join(parent_folder, patient)):
            niftis = glob.glob(os.path.join(parent_folder, patient, '*.nii.gz'))
            for nifti in niftis:
                file_path = os.path.split(nifti)[0]
                file_name = os.path.split(nifti)[1]
                if file_name[:len(patient)] != patient:
                    new_name = patient + '_' + file_name
                    os.rename(os.path.join(file_path, file_name), os.path.join(file_path, new_name))


def resample_CTs_in_patient_folders(nifti_path):
    contents = os.listdir(nifti_path)
    for patient in contents:
        niftis = glob.glob(os.path.join(nifti_path, patient, '*.nii.gz'))
        pet_path = ''
        ct_path = ''
        for nifti in niftis:
            filename = os.path.split(nifti)[1]
            if 'suv' in filename.lower():
                pet_path = nifti
            elif 'ctac' in filename.lower():
                ct_path = nifti
        if ct_path == '':
            print('No CTAC for {}'.format(patient))
        else:
            ct_new_path = ct_path[:-7] + '_resampled.nii.gz'
            resample_ct_nii_to_match_pet_nii(ct_path, pet_path, ct_new_path)




# dicom_path = os.path.join(r"\\onfnas01.uwhis.hosp.wisc.edu\radiology\Research\Bradshaw\Lymphoma_UW_Retrospective\Data\uw_analyzed", 'dicom' )
dicom_path = os.path.join(r"H:\Data\tmp", 'uw_analyzed' )

# nifti_path = os.path.join(r"\\onfnas01.uwhis.hosp.wisc.edu\radiology\Research\Bradshaw\Lymphoma_UW_Retrospective\Data\uw_analyzed", 'nifti')
# nifti_path = os.path.join(r"H:\Data\tmp", 'uw_analyzed_nifti')
nifti_path = os.path.join(r"\\onfnas01.uwhis.hosp.wisc.edu\radiology\Research\Bradshaw\Lymphoma_UW_Retrospective\Data\uw_analyzed", 'nifti')



#****************** One time things *********************
#convert the dataset
convert_all_files_to_nifti.convert_all_files_to_nifti(dicom_path, nifti_path, 'PT')

#resample the CTs to the PETs
resample_CTs_in_patient_folders(nifti_path)

#if you need to convert rtstructs for a specific annotator, use 'lee', 'mshin', or 'lee'
convert_specific_rtstruct_to_nifti(dicom_path, nifti_path, 'lee', 'PT')

#check which files were converted, save as csv
path_to_save_csv = os.path.join(nifti_path, 'converted_files.csv')
list_converted_nifti_files(nifti_path, path_to_save_csv)

#see if any subjects are missing a particular annotator, save in new csv
annotators = ['lee', 'mshin', 'scho']
for annotator in annotators:
    find_which_rtstructs_werent_converted(path_to_save_csv, annotator)
#****************** One time things *********************



#go through and create contours
annotators = ['lee', 'mshin', 'scho']
contents = os.listdir(nifti_path)
for annotator in annotators:
    #annotator-specific patients that have no lesions. need to create zero-mask files
    if annotator == 'lee':
        patients_with_no_lesions = ['4703', '4643', '4681', '4675', '4657', '4743', '4765', '4749', '4663', '4687', '4725'] #for Lee
    elif annotator == 'mshin':
        patients_with_no_lesions = ['4703', '4643', '4681', '4675', '4657', '4743', '4765', '4749', '4663']  #for MShin
    elif annotator == 'scho':
        patients_with_no_lesions = [ '4675', '4765'] #for Scho
    else:
        raise ValueError('No valid annotator specified')

    for patient in contents:
        print(patient)

        patient_dir = os.path.join(nifti_path, patient)
        if not os.path.isdir(patient_dir):
            continue

        all_folders  = get_all_terminal_subfolders(patient_dir)

        #if it's a patient with no contours, create empty mask
        no_lesions = 0
        for patient_no_lesion in patients_with_no_lesions:
            if patient_no_lesion in patient:
                no_lesions = 1
                all_files = os.listdir(patient_dir)
                pet = [s for s in all_files if 'SUV.nii' in s]
                roi_save_name1 = os.path.join(patient_dir, pet[0][:8] + '_no_lesions_combined' + annotator + '.nii.gz')
                roi_save_name2 = os.path.join(patient_dir, pet[0][:8] + '_no_lesions_combined_groups' + annotator + '.nii.gz')
                #create empty if it doesn't already exist
                if not os.path.isfile(roi_save_name1):
                    pet_nii = nib.load(os.path.join(patient_dir, pet[0]))
                    empty_mask = np.zeros(pet_nii.header.get_data_shape())
                    empty_nifti = nib.Nifti1Image(empty_mask, pet_nii.affine, pet_nii.header)
                    nib.save(empty_nifti, roi_save_name1 )
                    nib.save(empty_nifti, roi_save_name2 )
        #stop if patient has no lesions
        if no_lesions:
            continue

        #find the right roi folder
        roi_ind = []
        roi_folder = []
        for i, folder_i in enumerate(all_folders):
            if annotator in os.path.split(folder_i)[1].lower():
                roi_ind.append(i)
        if len(roi_ind) < 1:
            print("!!!****    No {} folders found in {}   ****!!!".format(annotator, patient))
        elif len(roi_ind) > 1:
            print("------ More than one roi folder found for {} in {} -------".format(annotator, patient))
            #too many mshins, one should say adj_mshin -- find it
            adj_ind = []
            for inds in roi_ind:
                if "adj" in os.path.split(all_folders[inds])[1].lower():
                    adj_ind.append(inds)

            if len(adj_ind) < 1:
                print("!!!****    Too many folders, no adj_ind found in {}   ****!!!".format(patient))
            elif len(adj_ind) > 1:
                print("!!!****    Too many adj folders found in {}   ****!!!".format(patient))
            else:
                roi_folder = all_folders[adj_ind[0]]
        else:
            roi_folder = all_folders[roi_ind[0]]


        ignore_strings = ['Struct_XD', 'Reference', 'Burden', 'DS1']
        group_strings = ['marrow', 'osseous', 'liver', 'extra-nodal', 'spleen', 'lymph-nodes']  #equivocal is already one
        try:
            combine_rois_into_one(roi_folder, patient_dir, patient, ignore_strings = ignore_strings, create_mip=True)
            combine_rois_into_groups(roi_folder, patient_dir, patient, group_strings=group_strings,
                                     ignore_strings = ignore_strings, create_mip=True)
        except:
            print('!!!!****  Unable to process this patient, {}  ****!!!!'.format(patient))







