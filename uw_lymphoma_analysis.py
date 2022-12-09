from convert_all_files_to_nifti import convert_all_files_to_nifti
import os
from get_all_terminal_subfolders import get_all_terminal_subfolders
from nifti_roi_tools import combine_rois_into_one, combine_rois_into_groups
import nibabel as nib
import numpy as np

# dicom_path = os.path.join(r"\\onfnas01.uwhis.hosp.wisc.edu\radiology\Research\Bradshaw\Lymphoma_UW_Retrospective\Data\uw_analyzed", 'dicom' )
dicom_path = os.path.join(r"H:\Data\tmp", 'uw_analyzed' )

# nifti_path = os.path.join(r"\\onfnas01.uwhis.hosp.wisc.edu\radiology\Research\Bradshaw\Lymphoma_UW_Retrospective\Data\uw_analyzed", 'nifti')
nifti_path = os.path.join(r"H:\Data\tmp", 'uw_analyzed_nifti')


#convert the dataset
convert_all_files_to_nifti.convert_all_files_to_nifti(dicom_path, nifti_path, 'PT')



#now create the combined ROI files

contents = os.listdir(nifti_path)
patients_with_no_lesions = ['4703', '4643', '4681', '4675', '4657', '4743', '4765', '4749', '4663']
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
            pet_nii = nib.load(os.path.join(patient_dir, pet[0]))
            empty_mask = np.zeros(pet_nii.header.get_data_shape())
            empty_nifti = nib.Nifti1Image(empty_mask, pet_nii.affine, pet_nii.header)
            nib.save(empty_nifti, os.path.join(patient_dir, pet[0][:8] + '_no_lesions_combined.nii.gz') )
            nib.save(empty_nifti, os.path.join(patient_dir, pet[0][:8] + '_no_lesions_combined_groups.nii.gz') )
    #stop if patient has no lesions
    if no_lesions:
        continue

    #find the right roi folder
    roi_ind = []
    roi_folder = []
    for i, folder_i in enumerate(all_folders):
        if "mshin" in os.path.split(folder_i)[1].lower():
            roi_ind.append(i)
    if len(roi_ind) < 1:
        print("!!!****    No mshin folders found in {}   ****!!!".format(patient))
    elif len(roi_ind) > 1:
        #too many mshins, one should say adj_mshin -- find it
        adj_ind = []
        for inds in roi_ind:
            if "adj" in os.path.split(all_folders[inds])[1].lower():
                adj_ind.append(inds)

        if len(adj_ind) < 1:
            print("!!!****    Too many mshin folders, no adj_ind found in {}   ****!!!".format(patient))
        elif len(adj_ind) > 1:
            print("!!!****    Too many adj_mshin folders found in {}   ****!!!".format(patient))
        else:
            roi_folder = all_folders[adj_ind[0]]
    else:
        roi_folder = all_folders[roi_ind[0]]


    ignore_strings = ['Struct_XD', 'Reference', 'Burden']
    group_strings = ['marrow', 'osseous', 'liver', 'extra-nodal', 'spleen', 'lymph-nodes']
    try:
        combine_rois_into_one(roi_folder, patient_dir, ignore_strings = ignore_strings, create_mip=True)
        combine_rois_into_groups(roi_folder, patient_dir, group_strings=group_strings, ignore_strings = ignore_strings, create_mip=True)
    except:
        print('!!!!****  Unable to process this patient, {}  ****!!!!'.format(patient))
