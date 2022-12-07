
import os
import glob
import nibabel as nib
from matplotlib import pyplot as plt
import numpy as np

rtstruct_folder = os.path.join(r"H:\Data\tmp\uw_analyzed_nifti\petlymph_4542_petlymph_4542\20150920_RTSTRUCT_ADJ_MShin")
save_folder = os.path.join(r"H:\Data\tmp\uw_analyzed_nifti\petlymph_4542_petlymph_4542")
ignore_strings = ['Struct_XD', 'Reference', 'Burden']
group_strings = ['marrow', 'osseous', 'liver', 'extra-nodal', 'spleen', 'eq', 'lymph-nodes']

def combine_rois_into_one(rtstruct_folder, save_folder, ignore_strings = None, create_mip=True):
    """
    rtstruct_folder - folder containing all nifti roi files
    save_folder - where to save combined mask
    ignore_strings - if filename contains one of these strings, it's ignored
    create_mip - create and save a png MIP of the combined mask

    """
    rt_niftis = glob.glob(os.path.join(rtstruct_folder, '*.nii.gz'))
    rt_combined = None

    for rt_i in rt_niftis:

        #ignore if includes ignore_strings
        ignore_rt = False
        for ignore in ignore_strings:
            if ignore.lower() in os.path.split(rt_i)[1].lower():
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


    #collapse values greater than 1
    rt_combined_data[rt_combined_data > 1] = 1
    combined_nifti = nib.Nifti1Image(rt_combined_data, rt_combined_affine, rt_combined_header)
    nib.save(combined_nifti, os.path.join(save_folder, os.path.split(rtstruct_folder)[1] + '_combined.nii.gz') )

    if create_mip:
        mip = np.max(rt_combined_data, axis=1)
        plt.imshow(mip)
        plt.imsave( os.path.join(save_folder, os.path.split(rtstruct_folder)[1] + '_combined_MIP.png'), mip)


def combine_rois_into_groups(rtstruct_folder, save_folder, group_strings, ignore_strings = None, create_mip=True):
    """
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
    dims = rt_initial.header.get_data_shape() + (len(group_strings),)

    roi_groups = np.zeros(dims)

    for rt_i in rt_niftis:
        group_ind = len(group_strings) - 1   #default is the nodal group
        #ignore if includes ignore_strings
        ignore_rt = False
        for ignore in ignore_strings:
            if ignore.lower() in os.path.split(rt_i)[1].lower():
                ignore_rt = True

        for i, groups in enumerate(group_strings):
            name_of_file = os.path.split(rt_i)[1].lower()
            if groups.lower() in name_of_file:
                if groups.lower() == 'eq':
                    if 'neq' in name_of_file or 'non-eq' in name_of_file:
                        #hmmm...does non-eq ever come before eq???? I'll ignore for now
                        group_ind = group_ind   #don't change
                    else:
                        #regardless of other label, add to equivocal group
                        group_ind = group_strings.index('eq')
                else:
                    group_ind = i

        #read image
        if not ignore_rt:
            rt_i_mask = nib.load(rt_i)
            rt_i_mask_data = rt_i_mask.get_fdata().copy()
            #add values to the correct group
            roi_groups[:,:,:,group_ind] = np.add(roi_groups[:,:,:,group_ind] , rt_i_mask_data)

    #collapse values greater than 1
    roi_groups[roi_groups > 1] = 1
    #go through each group, if it has data (value >1) create new nifti
    for i, groups in enumerate(group_strings):
        if np.max(roi_groups[:,:,:,i]) > 0:
            new_nifti_data = nib.Nifti1Image(roi_groups[:,:,:,i], rt_initial.affine, rt_initial.header)
            nib.save(new_nifti_data, os.path.join(save_folder, os.path.split(rtstruct_folder)[1] + '_combined_' + groups + '.nii.gz') )

            if create_mip:
                mip = np.max(new_nifti_data.get_fdata(), axis=1)
                plt.imshow(mip)
                plt.imsave( os.path.join(save_folder, os.path.split(rtstruct_folder)[1] + '_combined_' + groups + '_MIP.png'), mip)

