# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 19:57:48 2023

@author: PatrickBasitta
"""

import pandas as pd
import numpy as np
import csv
import os
import re
import argparse
import json

# using argparse for positinal arguments
parser = argparse.ArgumentParser()
parser.add_argument("-sf", "--qc_folder", type=str) 
parser.add_argument("-o", "--outdir", type=str)            
args = parser.parse_args()

directory = args.qc_folder

# get qc data
qc_files = os.listdir(directory)

# make dixt to get case idx
idx_qc_files_dict = {}
for idx, file in enumerate(qc_files):
    case_no = int(file.split("_")[0])
    idx_qc_files_dict[case_no] = idx

# sort case idx dict
sorted_idx_qc_files_dict = dict(sorted(idx_qc_files_dict.items()))
    
sorted_qc_files = []
for i in list(sorted_idx_qc_files_dict.values()):
    sorted_qc_files.append(qc_files[i])

# get qc data
qc_data_dict = {}
for index, qc_file in enumerate(sorted_qc_files):
    if qc_file.endswith(".csv"):
        if qc_file.split("_")[2] == "tumor":
            case_num = qc_file.split("_")[0]
        if qc_file.split("_")[2] == "normal":
            case_num = qc_file.split("_")[0]+"-N"
        QC_data = pd.read_csv(qc_file, delimiter= "\t", engine="python", header=None)
        QC_data = QC_data.T
        qc_data_dict[case_num] = list(QC_data.iloc[1])
        
# qc list
qc_data_lst = []

for i in range(len(qc_data_dict)):
    
    number = list(qc_data_dict.items())[i][0]
    lst_qc_data = list(qc_data_dict.items())[i][1]
    lst_qc_data.insert(0,number)
    qc_data_lst.append(lst_qc_data)
                            
# make dataframe 
header_col = QC_data.iloc[0]
header_col = list(header_col)
header_col.insert(0,"Sample")

# update qc_df
final_qc_data = pd.DataFrame(qc_data_lst, columns= header_col)  

# update qc_df
final_qc_data.insert(loc=0, column="Standort", value="Bonn")

# write csv
final_qc_data.to_csv(args.outdir+"qc.csv", index=False)