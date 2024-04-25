# -*- coding: utf-8 -*-
"""
Script to generate Patient DNPM data

Input:
    DatenmodellV2.0 xlsx of patient
    final processed variant data xlsx
    .
    .
    .

@author: Patrick Basitta
"""

import os
import pandas as pd
import json
import glob
import argparse
from collections import Counter
#import sys

# Using argparse for positinal arguments
parser = argparse.ArgumentParser(
    prog="DNPMdataset_generator")
parser.add_argument("-p", "--path", type=str)
args = parser.parse_args()

# global parameters
#PanCancer_to_PathoDB = "X:\PAT-Sequenzer\PANCANCER_OUTPUT_NXF\PanCancer_to_PathoDB"
PanCancer_to_PathoDB = args.path

# Get files an prefixes and add to lists
file_lst = []
prefix_lst = []
for file in os.listdir(PanCancer_to_PathoDB):
    if file.endswith(".xlsx"):  
        prefix_lst.append(file.split(".")[0])
        file_lst.append(file)
        
# check number of inputs files (2 files per patient)
check_num_prefixes = dict(Counter(prefix_lst))

for num_prefix in check_num_prefixes.values():
    if num_prefix != 2:
        raise ValueError(f"""2 files per patient expected! Only {str(num_prefix)} file in directory. 
            Make sure that all nessecary files are present!""")

# reduce duplicate prefixes       
prefix_lst = set(prefix_lst)

# process files
for prefix in prefix_lst:    
    file_dnpm_data_template = PanCancer_to_PathoDB+"/"+prefix+".meta.xlsx"

    # Get dnpm_data
    with open (file_dnpm_data_template, "rb") as f:
        dnpm_data = pd.read_excel(f,\
                             sheet_name=None, engine = "openpyxl")
        
    # get "Einfache Varianten"
    processed_variant_data = pd.read_excel(glob.glob(PanCancer_to_PathoDB+"/"+prefix+".*.final_processed.xlsx")[0])

    # get relevant data for report

    # get idx of reported
    idx_of_variants_to_report = processed_variant_data.index\
                             [processed_variant_data["submit"] == "x"].tolist()

    variants_to_report = processed_variant_data.loc[idx_of_variants_to_report, :]

    # get all necessary columns
    final_data = variants_to_report[["Chromosom",
                                    "Gen",
                                    "Transcript_ID",
                                    "Exon",
                                    "Position",
                                    "End Position",
                                    "Reference",
                                    "Allele",
                                    "cDNA_Change",
                                    "Amino_Acid_Change",
                                    "Coverage",
                                    "Frequency",
                                    "name dbsnp_v151_ensembl_hg38_no_alt_analysis_set",
                                     "Interpretation"]]

    # rename to DNPMdataset names
    final_data = final_data.rename(columns=\
               {"Reference": "Ref",
                "Allele": "Alt",
                "Coverage": "Read Depth",
                "Frequency": "Allelic Frequency", 
                "name dbsnp_v151_ensembl_hg38_no_alt_analysis_set": "dbSNP_ID",
                "CLIN_SIG": "Interpretation"})
        
    # add additional columns
    final_data.insert(loc=0, column="VariantenID", value="-")
    final_data.insert(loc=13, column="COSMIC ID", value="-")
                
    # VariantID-generator
    VariantID_lst = []
    for i in range(len(final_data)):
        i = i + 1
        VariantID_lst.append("VariantID_"+str(i))
        
    # insert VariantIDs
    final_data["VariantenID"] = VariantID_lst
    
    # add variant data to dnpm data
    dnpm_data["Einfache Variante"] = final_data
    
    #  make dictionary
    DNPM_tmp = dnpm_data["NGS-Bericht"].to_dict(orient='records')
    
    # convert dataframes into dictionaries
    dnpm_data_dict = {
        key: dnpm_data[key].to_dict(orient='records') 
        for key in dnpm_data.keys()
    }
    
    # add information
    DNPM_tmp[0]["Einfache Varianten"] = dnpm_data_dict["Einfache Variante"]
    DNPM_tmp[0]["Metadaten"] = dnpm_data_dict["NGS-Befund Metadaten"]
    DNPM_tmp[0]["Tumor Mutational Burden (TMB)"] = dnpm_data_dict\
                                                  ["Tumor Mutational Burden (TMB)"]
    
    def unlist_func(DNPM_tmp):
        
        for i in DNPM_tmp:
            return i
        
    DNPM_out = unlist_func(DNPM_tmp)
    
    # write to disk
    with open(prefix+".json", "w") as file:
        json.dump(
            DNPM_out, 
            file, 
            indent=4, 
            sort_keys=True
        )
    
    # Forward DNMP ID
    print(prefix)
    #print(prefix, file=sys.stdout)
    #sys.stdout.write(prefix)
    #sys.stdout.flush() 

# End    
###############################################################################
