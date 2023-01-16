
import os
import glob
import nibabel as nib
from nibabel.processing import resample_from_to
from matplotlib import pyplot as plt
import numpy as np
from scipy.ndimage import zoom, rotate


# rtstruct_folder = os.path.join(r"H:\Data\tmp\uw_analyzed_nifti\petlymph_4546_petlymph_4546\20150611_RTSTRUCT_ADJ_MShin")
# save_folder = os.path.join(r"H:\Data\tmp\uw_analyzed_nifti\petlymph_4546_petlymph_4546")
# ignore_strings = ['Struct_XD', 'Reference', 'Burden']
# group_strings = ['marrow', 'osseous', 'liver', 'extra-nodal', 'spleen', 'lymph-nodes']


def resample_and_normalize_nifti(nifti_img_path, new_voxel_dims, new_nifti_path, modality, crazy_rotation = False):
    #normalized to 0 mean, 1 std
    #modality can be PT, CT, or SUV
    #crazy rotation is to match the original deepmedic rotations that were used (crazy!)

    #threshold for measuring avg/std
    if modality == 'SUV':
        threshold = 0.1
    elif modality == 'PT':
        threshold = 500
    elif modality == 'CT':
        threshold = -300
    elif modality == 'RT':
        threshold = ''

    # Load the NIFTI image
    nii = nib.load(nifti_img_path)

    # Get the current voxel dimensions
    old_dims = nii.header.get_zooms()
    old_shape = nii.header.get_data_shape()

    # Calculate the scaling factors for the new voxel dimensions
    scaling_factors = (old_dims[0]/new_voxel_dims[0], old_dims[1]/new_voxel_dims[1], old_dims[2]/new_voxel_dims[2])

    #normalize image
    new_image = nii.get_fdata().copy()
    if threshold != '':
        avg = np.mean(new_image[new_image > threshold])
        std = np.std(new_image[new_image > threshold])
        new_image = (new_image - avg)/std


    # Resample the image
    if modality == 'RT':
        resampled_image = zoom(new_image, scaling_factors, order=0) #order 0 = NN ...i think
    else:
        resampled_image = zoom(new_image, scaling_factors, order=1) #order 1 = linear ... i think
    new_shape = resampled_image.shape



    # Create a new NIFTI image with the resampled data and the updated voxel dimensions
    # !! need to fix this affine -- DOESN'T CONSIDER DIRECTION, AND assumes the image is centered on 0, which isn;t
    # always the case!!
    new_start_x = new_voxel_dims[0]*new_shape[0]/2
    new_start_y = new_voxel_dims[1]*new_shape[1]/2
    #z isn't usually centered at 0
    old_z_start = nii.affine[2,3]
    if nii.affine[2,2] > 0:  #defines direction which to add extra
        new_start_z = old_z_start - ( new_voxel_dims[2]*new_shape[2] - old_dims[2]*old_shape[2] )/2
    else:
        new_start_z = old_z_start + ( new_voxel_dims[2]*new_shape[2] - old_dims[2]*old_shape[2] )/2

    #new affine
    nii_affine_new = nii.affine.copy()
    nii_affine_new[0,0] = np.sign(nii.affine[0,0]) * new_voxel_dims[0]
    nii_affine_new[1,1] = np.sign(nii.affine[1,1]) * new_voxel_dims[1]
    nii_affine_new[2,2] = np.sign(nii.affine[2,2]) * new_voxel_dims[2]
    nii_affine_new[:3,3] = [new_start_x, new_start_y, new_start_z]
    #new header
    nii_header_new = nii.header.copy()
    nii_header_new['pixdim'][1:4] = new_voxel_dims

    if crazy_rotation:
        resampled_image = rotate(resampled_image, 90)

    resampled_nii = nib.Nifti1Image(resampled_image, nii_affine_new, nii_header_new)

    # Save the resampled image
    nib.save(resampled_nii, new_nifti_path)


def resample_nifti_to_pixel_dims(nifti_path, new_voxel_dims, new_nifti_path):
    # Load the NIFTI image
    nii = nib.load(nifti_path)

    # Get the current voxel dimensions
    old_dims = nii.header.get_zooms()
    old_shape = nii.header.get_data_shape()

    # Calculate the scaling factors for the new voxel dimensions
    scaling_factors = (old_dims[0]/new_voxel_dims[0], old_dims[1]/new_voxel_dims[1], old_dims[2]/new_voxel_dims[2])

    # Resample the image
    resampled_image = zoom(nii.get_fdata(), scaling_factors, order=1)
    new_shape = resampled_image.shape

    # Create a new NIFTI image with the resampled data and the updated voxel dimensions
    # !! need to fix this affine -- DOESN'T CONSIDER DIRECTION, AND assumes the image is centered on 0, which isn;t
    # always the case!!
    new_start_x = new_voxel_dims[0]*new_shape[0]/2
    new_start_y = new_voxel_dims[1]*new_shape[1]/2
    #z isn't usually centered at 0
    old_z_start = nii.affine[2,3]
    if nii.affine[2,2] > 0:  #defines direction which to add extra
        new_start_z = old_z_start - ( new_voxel_dims[2]*new_shape[2] - old_dims[2]*old_shape[2] )/2
    else:
        new_start_z = old_z_start + ( new_voxel_dims[2]*new_shape[2] - old_dims[2]*old_shape[2] )/2

    #new affine
    nii_affine_new = nii.affine.copy()
    nii_affine_new[0,0] = np.sign(nii.affine[0,0]) * new_voxel_dims[0]
    nii_affine_new[1,1] = np.sign(nii.affine[1,1]) * new_voxel_dims[1]
    nii_affine_new[2,2] = np.sign(nii.affine[2,2]) * new_voxel_dims[2]
    nii_affine_new[:3,3] = [new_start_x, new_start_y, new_start_z]
    #new header
    nii_header_new = nii.header.copy()
    nii_header_new['pixdim'][1:4] = new_voxel_dims
    resampled_nii = nib.Nifti1Image(resampled_image, nii_affine_new, nii_header_new)

    # Save the resampled image
    nib.save(resampled_nii, new_nifti_path)


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
    lesion (e.g., 1=marrow, 2=osseous, ...) For equivocal (EQ) lesions, these are assigned a number + total number of
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
                group_ind = i+1  #plus 1, or else first one will be zero!!!

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






