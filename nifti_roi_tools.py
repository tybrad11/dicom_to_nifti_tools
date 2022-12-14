
import os
import glob
import nibabel as nib
from nibabel.processing import resample_from_to
from matplotlib import pyplot as plt
import numpy as np

# rtstruct_folder = os.path.join(r"H:\Data\tmp\uw_analyzed_nifti\petlymph_4546_petlymph_4546\20150611_RTSTRUCT_ADJ_MShin")
# save_folder = os.path.join(r"H:\Data\tmp\uw_analyzed_nifti\petlymph_4546_petlymph_4546")
# ignore_strings = ['Struct_XD', 'Reference', 'Burden']
# group_strings = ['marrow', 'osseous', 'liver', 'extra-nodal', 'spleen', 'lymph-nodes']

def resample_nifti_to_pixel_dims(nifti_path, voxel_dims)



def resample_ct_nii_to_match_pet_nii(ct_path, pt_path, ct_save_path):
    ct_i = nib.load(ct_path)
    pt_i = nib.load(pt_path)

    #resample
    ct_rs_i = resample_from_to(ct_i, pt_i)

    #saving
    nib.save(ct_rs_i, ct_save_path)

def determine_if_rt_should_be_included(name_of_file):
    #yd means scho deletes, xd means mshin deletes.
    name_of_file = name_of_file.replace(',','_').lower()
    if 'yd_' in name_of_file:
        return False
    if 'xd_' in name_of_file and 'yl_' not in name_of_file:  #yl overrides xd
        return False
    return True

def determine_if_rt_is_equivocal(name_of_file):
    name_of_file = name_of_file.replace(',','_').lower()
    if 'neq' in name_of_file or 'non-eq' in name_of_file:
        if '_eq' not in name_of_file:
            return False
        else:
            #neq or non-eq might come before or after eq
            ind_neq = name_of_file.find('neq')
            ind_noneq = name_of_file.find('non_eq')
            ind_eq = name_of_file.find('_eq')
            #neq after eq means neq
            if max(ind_noneq, ind_neq) > ind_eq:
                return False
            else: #neq before eq means eq
                return True

    elif '_eq' in name_of_file:
        return True
    else:
        return False



def combine_rois_into_one(rtstruct_folder, save_folder, patient, ignore_strings = None, create_mip=True):
    """
    Combine all rois in the separate nifti files into a single binary mask, ignoring filenames that contain certain
    strings. Save as a nifti and MIP.png
    rtstruct_folder - folder containing all nifti roi files
    save_folder - where to save combined mask
    ignore_strings - if filename contains one of these strings, it's ignored
    annotator - string indicating which annotator it was, added to filenames
    create_mip - create and save a png MIP of the combined mask

    """
    rt_niftis = glob.glob(os.path.join(rtstruct_folder, '*.nii.gz'))
    rt_combined = None

    for rt_i in rt_niftis:

        #ignore if includes ignore_strings
        ignore_rt = False
        name_of_file = os.path.split(rt_i)[1].lower()
        for ignore in ignore_strings:
            if ignore.lower() in name_of_file:
                ignore_rt = True

        if not determine_if_rt_should_be_included(name_of_file):
            ignore_rt = True

        #read image
        if not ignore_rt:
            rt_i_mask = nib.load(rt_i)

            #if it's the first one, set as combined mask
            if rt_combined == None:
                rt_combined  = 1
                rt_combined_data = rt_i_mask.get_fdata().copy()
                rt_combined_header = rt_i_mask.header
                rt_combined_affine = rt_i_mask.affine
            else:
                rt_combined_data = np.add(rt_combined_data,  rt_i_mask.get_fdata().copy())
        else:
            continue


    #collapse values greater than 1
    rt_combined_data[rt_combined_data > 1] = 1
    combined_nifti = nib.Nifti1Image(rt_combined_data, rt_combined_affine, rt_combined_header)
    nib.save(combined_nifti, os.path.join(save_folder, patient + '_' + os.path.split(rtstruct_folder)[1] + '_combined.nii.gz') )

    if create_mip:
        mip = np.max(rt_combined_data, axis=1)
        plt.imshow(mip)
        plt.imsave( os.path.join(save_folder, patient + '_' + os.path.split(rtstruct_folder)[1] + '_combined_MIP.png'), mip)




def combine_rois_into_groups(rtstruct_folder, save_folder, patient, group_strings, ignore_strings = None, create_mip=True):
    """
    Combine all rois in the separate nifti files into a single mask where each number represents a different class of
    lesion (e.g., 0=marrow, 1=osseous, ...) For equivocal (EQ) lesions, these are assigned a number + total number of
    groups. E.g., if liver = 2, and number of groups = 7, liver EQ is 2 + 7. Also ignore files that contain strings in
    ignore_strings.

    rtstruct_folder - folder containing all nifti roi files
    save_folder - where to save combined mask
    group_strings - which groups of ROIs to create
    ignore_strings - if filename contains one of these strings, it's ignored
    create_mip - create and save a png MIP of the combined mask

    """
    rt_niftis = glob.glob(os.path.join(rtstruct_folder, '*.nii.gz'))

    #initialize array of groups, one for each group plus a "regular"
    #read in first dataset to get header info
    rt_initial = nib.load(rt_niftis[0])
    dims = rt_initial.header.get_data_shape()

    roi_groups = np.zeros(dims)

    for rt_i in rt_niftis:
        name_of_file = os.path.split(rt_i)[1].lower()
        group_ind = len(group_strings) - 1   #default is the nodal group
        #ignore if includes ignore_strings
        ignore_rt = False
        for ignore in ignore_strings:
            if ignore.lower() in os.path.split(rt_i)[1].lower():
                ignore_rt = True

        if not determine_if_rt_should_be_included(name_of_file):
            ignore_rt = True

        #see what group it belongs to based on filename
        for i, groups in enumerate(group_strings):
            if groups.lower() in name_of_file:
                group_ind = i

        #see if it's equivocal, if it is add as it's own category (each group has it owns eq category)
        is_equivocal = determine_if_rt_is_equivocal(name_of_file)
        if is_equivocal:
            group_ind = group_ind + len(group_strings)

        #read image, set mask to correct group number
        if not ignore_rt:
            rt_i_mask = nib.load(rt_i)
            rt_i_mask_data = rt_i_mask.get_fdata().copy()
            #add values to the correct group
            roi_groups[rt_i_mask_data > 0] = group_ind


    new_nifti_data = nib.Nifti1Image(roi_groups, rt_initial.affine, rt_initial.header)
    nib.save(new_nifti_data, os.path.join(save_folder, patient + '_' + os.path.split(rtstruct_folder)[1] + '_combined_groups.nii.gz') )

    if create_mip:
        mip = np.max(new_nifti_data.get_fdata(), axis=1)
        plt.imshow(mip)
        plt.imsave( os.path.join(save_folder, patient + '_' + os.path.split(rtstruct_folder)[1] + '_combined_groups_MIP.png'), mip)






