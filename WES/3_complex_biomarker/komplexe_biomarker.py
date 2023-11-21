# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 14:05:54 2023

@author: 50118212
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
parser.add_argument("-hf", "--hrd_folder", type=str) 
parser.add_argument("-cp", "--cp_folder", type=str)
parser.add_argument("-o", "--outdir", type=str)            
args = parser.parse_args()

# make dataframe !!!check Columns before submittung!!!
df_komplexe_b = pd.DataFrame()

header_col = ["Standort", "Sample", "Anzahl TMB Mut. missense", \
              "Anzahl TMB Mut. InDel", "Regionsgroesse [Mb]", "TMB Missense", \
              "TMB Missense + InDel", "LOH", "TAI", "LST", "HRD", "TZG", \
              "Ploidie", "MSI Anzahl sites", "MSI Anzahl instabile sites", \
              "MSI"]

for col in header_col:
    df_komplexe_b[col] = ""

# HRD
for index, clc_file_hrd in enumerate(os.listdir\
            (args.hrd_folder)):
    if clc_file_hrd.endswith(".json"):
        
        # get case number
        if clc_file_hrd.split("_")[2] != "1":
             case_num = clc_file_hrd.split("_")[2]
             df_komplexe_b.loc[index,"Sample"] = int(case_num)
             
        # Standort
        df_komplexe_b.loc[index,"Standort"] = "Bonn"
                   
        # process files
        with open(args.hrd_folder+clc_file_hrd) as hrd_file:
            HRD_data = hrd_file.read()
        #print(HRD_data)
        parsed_HRD_data = json.loads(HRD_data)
        #print(parsed_HRD_data["data"]["hrd_score"]["table_1"][0]["id"])
        #print(parsed_HRD_data["data"]["hrd_score"]["table_1"][0]["value"])

        for i in range(len(parsed_HRD_data["data"]["hrd_score"]["table_1"])):
            if parsed_HRD_data["data"]["hrd_score"]["table_1"][i]["id"] == "HRD score":
                df_komplexe_b.loc[index,"HRD"] = parsed_HRD_data\
                                        ["data"]["hrd_score"]["table_1"][i]["value"]
            if parsed_HRD_data["data"]["hrd_score"]["table_1"][i]["id"] == "LOH score":
                df_komplexe_b.loc[index,"LOH"] = parsed_HRD_data\
                                        ["data"]["hrd_score"]["table_1"][i]["value"]
            if parsed_HRD_data["data"]["hrd_score"]["table_1"][i]["id"] == "LST score":
               df_komplexe_b.loc[index,"LST"] = parsed_HRD_data\
                                        ["data"]["hrd_score"]["table_1"][i]["value"]
            if parsed_HRD_data["data"]["hrd_score"]["table_1"][i]["id"] == "TAI score":
                df_komplexe_b.loc[index,"TAI"] = parsed_HRD_data\
                                        ["data"]["hrd_score"]["table_1"][i]["value"]

# Celluarity and Ploidy
for index, seqz_file_CP in enumerate(os.listdir\
            (args.cp_folder)):
    if seqz_file_CP.endswith(".txt"):
       
        # get case number
        if seqz_file_CP.split("_")[0] != "19":
             case_num_CP = seqz_file_CP.split("_")[0]
        
        # process files
        with open(args.cp_folder+seqz_file_CP) as CP_file:
            CP_data = pd.read_csv(CP_file, delimiter= "\t", engine="python")
            celluarity_data = round(CP_data["cellularity"].iloc[0]*100,2)
            ploidy_data = CP_data["ploidy"].iloc[0]
        
        index_CP = df_komplexe_b.index[df_komplexe_b["Sample"]==int(case_num_CP)]
        df_komplexe_b["TZG"].iloc[index_CP] = celluarity_data
        df_komplexe_b["Ploidie"].iloc[index_CP] = ploidy_data 
  
# MSI
# TMB

# Final output
df_komplexe_b = df_komplexe_b.sort_values(by=["Sample"], ascending=True)    
df_komplexe_b.to_csv(args.outdir+"biomarker.csv", index=False)






