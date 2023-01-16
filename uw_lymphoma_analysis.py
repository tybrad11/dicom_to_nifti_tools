from convert_all_files_to_nifti import convert_all_files_to_nifti, convert_specific_rtstruct_to_nifti
import os
from get_all_terminal_subfolders import get_all_terminal_subfolders
from nifti_roi_tools import combine_rois_into_one, combine_rois_into_groups, resample_ct_nii_to_match_pet_nii, resample_and_normalize_nifti
import nibabel as nib
import numpy as np
from list_converted_nifti_files import list_converted_nifti_files, find_which_rtstructs_werent_converted
import glob

#NOTE: evaluation of results after deepemedic is performed in Matlab (Amy Weisman's scripts): C:\Users\tjb129\Box Sync\MATLAB\lymphoma_uw_retrospective

#NOTE: group_strings = ['marrow', 'osseous', 'liver', 'extra-nodal', 'spleen', 'lymph-nodes']
#meaning 1='marrow', 2= 'osseous',3= 'liver',4= 'extra-nodal',5= 'spleen',6= 'lymph-nodes'
#actually, I think this got messed up and starts at 0=Marrow...ugh. Redo

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





## go through and create contours for each annotator ##

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
        print(patient + ', ' + annotator)

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
                roi_save_name1 = os.path.join(patient_dir, patient + '_no_lesions_combined_' + annotator + '.nii.gz')
                roi_save_name2 = os.path.join(patient_dir, patient + '_no_lesions_combined_groups_' + annotator + '.nii.gz')
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
            combine_rois_into_one(roi_folder, patient_dir, patient, ignore_strings, create_mip=True)
            combine_rois_into_groups(roi_folder, patient_dir, patient, group_strings=group_strings,
                                     ignore_strings = ignore_strings, create_mip=True)
        except:
            print('!!!!****  Unable to process this patient, {}  ****!!!!'.format(patient))



## precprocess for deepmedic - resample everything to 2x2x2 and normalize ##
deepmedic_save_path = os.path.join(r"\\onfnas01.uwhis.hosp.wisc.edu\radiology\Research\Bradshaw\Lymphoma_UW_Retrospective\Data\uw_analyzed", 'for_deepmedic')
nifti_path = os.path.join(r"\\onfnas01.uwhis.hosp.wisc.edu\radiology\Research\Bradshaw\Lymphoma_UW_Retrospective\Data\uw_analyzed", 'nifti')

annotators = ['lee', 'mshin', 'scho']
contents = os.listdir(nifti_path)
rs_vox_dims = [2,2,2] # in mm

for patient in contents:
    print(patient)
    niftis = glob.glob(os.path.join(nifti_path, patient, '*.nii.gz'))
    pet_niftis = [i for i in niftis if 'SUV.nii' in i]
    ct_niftis = [i for i in niftis if ('CTAC' in i and 'resampled' not in i)]
    if len(pet_niftis) < 1:
        print("Too few PETs for {}".format(patient))
        continue
    if len(ct_niftis) < 1:
        print("Too few CTs for {}".format(patient))
        continue

    new_pet_path = os.path.join(deepmedic_save_path, patient + '_PET_norm_' + str(rs_vox_dims[0]) + 'x' + str(rs_vox_dims[1]) + 'x' + str(rs_vox_dims[0]) + '_rotated.nii.gz')
    new_ct_path = os.path.join(deepmedic_save_path, patient + '_CTAC_norm_' + str(rs_vox_dims[0]) + 'x' + str(rs_vox_dims[1]) + 'x' + str(rs_vox_dims[0]) + '_rotated.nii.gz')

    resample_and_normalize_nifti(pet_niftis[0], rs_vox_dims, new_pet_path, 'SUV', crazy_rotation = True)
    resample_and_normalize_nifti(ct_niftis[0], rs_vox_dims, new_ct_path, 'CT', crazy_rotation = True)


    for annotator in annotators:
        roi_niftis = [i for i in niftis if (annotator in i.lower() and 'combined_groups' in i)]
        if len(roi_niftis) < 1:
            print('Too few ROIs for subject/annotator {} , {} -- selecting first one:'.format(patient, annotator))
            continue

        new_roi_path = os.path.join(deepmedic_save_path, patient + '_RTSTRUCT_' + annotator + '_groups_' + str(rs_vox_dims[0]) + 'x' + str(rs_vox_dims[1]) + 'x' + str(rs_vox_dims[0]) + '_rotated.nii.gz')
        resample_and_normalize_nifti(roi_niftis[0], rs_vox_dims, new_roi_path, 'RT', crazy_rotation = True)


# create test list for deepmedic
annotators = ['lee', 'mshin', 'scho']
template_path = '/Data/Bradshaw/Lymphoma_UW_Retrospective/Data/uw_analyzed/for_deepmedic/'
save_file_path = os.path.join(r"\\onfnas01.uwhis.hosp.wisc.edu\radiology\Research\Bradshaw\Lymphoma_UW_Retrospective\Analysis\deepmedic_files")
deepmedic_path = os.path.join(r"\\onfnas01.uwhis.hosp.wisc.edu\radiology\Research\Bradshaw\Lymphoma_UW_Retrospective\Data\uw_analyzed", 'for_deepmedic')
patient_string_length = len('petlymph_1234_petlymph_1234')
save_pt_file = 'testChannels_PET.cfg'
save_ct_file = 'testChannels_CT.cfg'
save_label_file = 'testGtLabels_{}.cfg'
save_name_file = 'testNamesOfPredictions.cfg'

contents = os.listdir(deepmedic_path)
contents = sorted(contents)
patients = []

for folders in contents:
    patient = folders[:patient_string_length]
    #only include new patients and even-numbered patients (baseline scans)
    if int(patient[-1]) % 2 == 0:
        if patient not in patients:
            patients.append(patient)
            pets =  [i for i in contents if ('_PET_' in i and patient in i)]
            cts =  [i for i in contents if ('_CTAC_' in i and patient in i)]
            with open(os.path.join(save_file_path, save_pt_file), 'a') as f:
                f.write(template_path + pets[0] + '\n')
            with open(os.path.join(save_file_path, save_ct_file), 'a') as f:
                f.write(template_path + cts[0] + '\n')
            with open(os.path.join(save_file_path, save_name_file), 'a') as f:
                f.write(template_path + patient + '\n')
            for annotator in annotators:
                label = [i for i in contents if ('RTSTRUCT_' + annotator in i and patient in i)]
                with open(os.path.join(save_file_path, save_label_file.format(annotator)), 'a') as f:
                    f.write(template_path + label[0] + '\n')




# combine the ensemble of segmentation masks
parent_folder = os.path.join(r"\\onfnas01.uwhis.hosp.wisc.edu\radiology\Research\Bradshaw\Lymphoma_UW_Retrospective\Data\deepmedic_model_predictions")
segmentation_folder1 = os.path.join(parent_folder, 'ensemble1_MSSM_WB_Fold1')
segmentation_folder2 = os.path.join(parent_folder, 'ensemble2_MSSM_WB_Fold2')
segmentation_folder3 = os.path.join(parent_folder, 'ensemble3_MSSM_WB_Fold3')
# segmentation_folder1 = os.path.join(parent_folder, 'ensemble1_updated')
# segmentation_folder2 = os.path.join(parent_folder, 'ensemble2_updated')
# segmentation_folder3 = os.path.join(parent_folder, 'ensemble3_updated')
output_folder = os.path.join(parent_folder, 'intersection_segmentation_maps')

contents = os.listdir(segmentation_folder1)
contents = sorted(contents)

for patient in contents:
    if patient[-7:] != '.nii.gz':
        continue
    s1 = nib.load(os.path.join(segmentation_folder1, patient))
    s2 = nib.load(os.path.join(segmentation_folder2, patient))
    s3 = nib.load(os.path.join(segmentation_folder3, patient))

    s_sum = s1.get_fdata() + s2.get_fdata() + s3.get_fdata()

    s_sum[s_sum < 3] = 0
    s_sum[s_sum > 2] = 1

    s_final = nib.Nifti1Image(s_sum, s1.affine, s1.header)
    nib.save(s_final, os.path.join(output_folder, patient))



# fix problemeatic ones -- list the dimensions to see how to fix
pts= ['4672','4744','4754','4780'  ,  '4718'  ,  '4676','4688', '4710']
for pt_num in pts:
    patient = 'petlymph_{}_petlymph_{}'.format(pt_num, pt_num)
    print(patient)
    orig_pat_folder = os.path.join(r"\\onfnas01.uwhis.hosp.wisc.edu\radiology\Research\Bradshaw\Lymphoma_UW_Retrospective\Data\uw_analyzed\nifti", patient)

    input_pt = nib.load(os.path.join(deepmedic_path, patient + '_PET_norm_2x2x2_rotated.nii.gz'))
    input_ct = nib.load(os.path.join(deepmedic_path, patient + '_CTAC_norm_2x2x2_rotated.nii.gz'))
    dm_seg = nib.load(os.path.join(output_folder, patient + '_Segm.nii.gz'))
    lee = nib.load(os.path.join(deepmedic_path, patient + '_RTSTRUCT_lee_groups_2x2x2_rotated.nii.gz'))
    mshin = nib.load(os.path.join(deepmedic_path, patient + '_RTSTRUCT_mshin_groups_2x2x2_rotated.nii.gz'))
    scho = nib.load(os.path.join(deepmedic_path, patient + '_RTSTRUCT_scho_groups_2x2x2_rotated.nii.gz'))

    niftis = glob.glob(os.path.join(orig_pat_folder, '*.nii.gz'))
    for n in niftis:
        n_name = os.path.split(n)[1]
        if 'SUV.nii.gz' in n_name:
            orig_pt_path= n
        if 'CTAC.nii.gz' in n_name:
            orig_ct_path= n
        elif 'lee' in n_name or 'Lee' in n_name or 'LEE' in n_name:
            orig_lee_path = n
        elif 'scho' in n_name or 'SCho' in n_name or 'SCHO' in n_name:
            orig_scho_path = n

    orig_pt = nib.load(orig_pt_path)
    orig_ct = nib.load(orig_ct_path)
    orig_lee =  nib.load(orig_lee_path)
    orig_scho =  nib.load(orig_scho_path)
    print('processed_pet_and_ct_dm')
    print(input_pt.header.get_data_shape())
    print(input_ct.header.get_data_shape())
    print('dm_output')
    print(dm_seg.header.get_data_shape())
    print('processed_masks_dm')
    print(lee.header.get_data_shape())
    print(mshin.header.get_data_shape())
    print(scho.header.get_data_shape())

    print('origs')
    print(orig_pt.header.get_data_shape())
    print(orig_ct.header.get_data_shape())
    print(orig_lee.header.get_data_shape())
    print(orig_scho.header.get_data_shape())
    print('--')

#
# # delete the extra (old) combined .nii now that we have new ones. 'lee.nii.gz', 'mshin.nii.gz'
# nifti_path = os.path.join(r"\\onfnas01.uwhis.hosp.wisc.edu\radiology\Research\Bradshaw\Lymphoma_UW_Retrospective\Data\uw_analyzed", 'nifti')
#
# annotators = ['lee', 'mshin', 'scho']
# contents = os.listdir(nifti_path)
#
# for patient in contents:
#     print(patient)
#     niftis = glob.glob(os.path.join(nifti_path, patient, '*.nii.gz'))
#     for annotator in annotators:
#         roi_niftis = [i for i in niftis if (annotator in i.lower() and 'combined_groups' in i)]
#
#         if len(roi_niftis) > 1:
#             bad_ending = annotator + '.nii.gz'
#             for r in roi_niftis:
#
#                 if bad_ending in r:
#                     print(r)
#                     os.remove(r)
#                     print('Deleted {}'.format(r))
#

