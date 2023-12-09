# -*- coding: utf-8 -*-
"""
Created on Tue Nov 21 10:57:32 2023
@author: Patrick Basitta
Script for collecting CNV data in accordance to the WES Pilot specification
"""

import pandas as pd
import os
import argparse

# using argparse for positinal arguments
parser = argparse.ArgumentParser()
parser.add_argument("-hf", "--cnv_folder", type=str)  
parser.add_argument("-o", "--outdir", type=str)         
args = parser.parse_args()

# make dataframe
df_cnv = pd.DataFrame()

header_col = ["Standort", "Sample", "Chr", "Start", "Ende", "CN", \
              "CN_A", "CN_B"]

for col in header_col:
    df_cnv[col] = ""

# sort file list after numbers    
folder_lst = os.listdir(args.cnv_folder)
order_dict = {}
for idx, file in enumerate(folder_lst):
    sort_num = file.split("_")[0]
    order_dict[idx] = int(sort_num)

sorted_order_dict = sorted(order_dict.items(), key=lambda x:x[1])

idx_lst = []
for case_idx in sorted_order_dict:
    idx_lst.append(case_idx[0])

orderd_file_lst = [folder_lst[i] for i in idx_lst]
    

for index, cnv_segments in enumerate(orderd_file_lst):
    if cnv_segments.endswith(".txt"):
        
        # get case number
        case_num_cnv = int(cnv_segments.split("_")[0])
               
        # process files
        with open(args.cnv_folder+cnv_segments) as cnv_file:
            cnv_data = pd.read_csv(cnv_file, delimiter= "\t", engine="python")
            
            # rename columns
            cnv_data = cnv_data.rename(columns={"chromosome": "Chr",
                                               "start.pos": "Start",
                                               "end.pos": "Ende",
                                               "CNt": "CN",
                                               "A": "CN_A",
                                               "B": "CN_B"})
            
            # get case specif columns
            cnv_data_case = cnv_data[["Chr", "Start", "Ende", "CN", "CN_A", "CN_B"]]
            
            # Standort
            cnv_data_case.insert(loc=0, column="Standort", value="Bonn") 
            cnv_data_case.insert(loc=1, column="Sample", value=case_num_cnv)
            
            dfs = [df_cnv, cnv_data_case]
            df_cnv = pd.concat(dfs, ignore_index=True)

df_cnv["Chr"] = "chr" + df_cnv["Chr"].astype(str)
            
df_cnv.to_csv(args.outdir+"cnv.csv", index=False)
            
            
            
            
          
