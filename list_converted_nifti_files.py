
import os
import pandas as pd

parent_folder_after_nifti_conversion = os.path.join(r"H:\Data\tmp\uw_analyzed_nifti")
path_to_save_csv = os.path.join(r"H:\Data\tmp\uw_analyzed_nifti\converted_files.csv")

output = pd.DataFrame(columns = ['Subject', 'PT', 'PT_SUV', 'CT', 'RTSTRUCT-Name', 'RTSTRUCT-ROI'])

for subj_folder in os.listdir(parent_folder_after_nifti_conversion):
    subj_path = os.path.join(parent_folder_after_nifti_conversion, subj_folder)
    subfolders = []
    niftis = []

    #store all nifti files and subfolders into lists
    for files_or_folder in os.listdir(subj_path):
        file_path = os.path.join(subj_path, files_or_folder)
        if files_or_folder.endswith('.nii.gz'):
            niftis.append(files_or_folder)
        elif os.path.isdir(file_path):
            subfolders.append(file_path)
    #see if subfolders have subfolders.
    sub_subfolders = []
    for subfolder_i in subfolders:
        subfolder_path = os.path.join(subj_path, subfolder_i)
        for file_i in os.listdir(subfolder_path):
            file_path = os.path.join(subfolder_path, file_i)
            if os.path.isdir(file_path):
                sub_subfolders.append(file_path)
    subfolders = subfolders + sub_subfolders

    #get all niftis in main folder, check for PET and CT
    for nifti_i in niftis:
        if 'SUV.nii' in nifti_i:
            output = output.append({'Subject': subj_folder,
                                    'PT': '-',
                                    'PT_SUV': nifti_i,
                                    'CT' : '-',
                                    'RTSTRUCT-Name': '-',
                                    'RTSTRUCT-ROI' : '-'}, ignore_index=True)
        elif 'PT_' in nifti_i or 'MAC' in nifti_i:
            output = output.append({'Subject': subj_folder,
                                    'PT': nifti_i,
                                    'PT_SUV': '-',
                                    'CT' : '-',
                                    'RTSTRUCT-Name': '-',
                                    'RTSTRUCT-ROI' : '-'}, ignore_index=True)
        elif 'CT_' in nifti_i or 'CTAC' in nifti_i:
            output = output.append({'Subject': subj_folder,
                                    'PT': '-',
                                    'PT_SUV': '-',
                                    'CT' : nifti_i,
                                    'RTSTRUCT-Name': '-',
                                    'RTSTRUCT-ROI' : '-'}, ignore_index=True)

    #now get all the RTSTRUCT ROIS in the subfoldres
    for subfolder_i in subfolders:
        rt_niftis = []
        subfolder_path = os.path.join(subj_path, subfolder_i)
        for file_i in os.listdir(subfolder_path):
            file_path = os.path.join(subfolder_path, file_i)
            if file_i.endswith('.nii.gz'):
                rt_niftis.append(file_i)
        if len(rt_niftis) > 0:
            for roi_i in rt_niftis:
                output = output.append({'Subject': subj_folder,
                                    'PT': '-',
                                    'PT_SUV': '-',
                                    'CT' : '-',
                                    'RTSTRUCT-Name': subfolder_i,
                                    'RTSTRUCT-ROI' : roi_i}, ignore_index=True)

output.to_csv(path_to_save_csv)
