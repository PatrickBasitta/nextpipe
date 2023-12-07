# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 08:05:16 2023

@author: PatrickBasitta
"""

import pandas as pd
import numpy as np
import csv
import os
import re
import argparse
import json
# import glob

# using argparse for positinal arguments
parser = argparse.ArgumentParser()
parser.add_argument("-sf", "--snv_folder", type=str) 
parser.add_argument("-o", "--outdir", type=str)            
args = parser.parse_args()

# get csv files list from a folder
snv_folder = args.snv_folder
#snv_folder = "WES_Pilot_final/somatic_variants/csv"
# snv_files = glob.glob(snv_folder+ "/*.csv")
snv_files = os.listdir(snv_folder)

# make dixt to get case idx
idx_snv_files_dict = {}
for idx, file in enumerate(snv_files):
    case_no = int(file.split(".")[0])
    idx_snv_files_dict[case_no] = idx

# sort case idx dict
sorted_idx_snv_files_dict = dict(sorted(idx_snv_files_dict.items()))
    
sorted_snv_files = []
for i in list(sorted_idx_snv_files_dict.values()):
    sorted_snv_files.append(snv_files[i])
    
    
# read each file to a df
tmp_dfs_lst = (pd.read_csv(snv_file) for snv_file in sorted_snv_files)

# concat all dfs
final_df = pd.concat(tmp_dfs_lst, ignore_index=True)

# add chr prefix
final_df["Chr"] = "chr" + final_df["Chr"].astype(str)

# save to csv
final_df.to_csv(args.outdir+"varianten.csv", index=False)


