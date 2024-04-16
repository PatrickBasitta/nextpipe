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
import random


# Get template file
file_dnpm_data_template="PATIENT_DNPM_BD-0001-DatenmodellV2.0.xlsx"

# Get dnpm_data
with open (file_dnpm_data_template, "rb") as f:
    dnpm_data = pd.read_excel(f,\
                         sheet_name=None, engine = "openpyxl")
        
        
# get "Einfache Varianten"
processed_variant_data = pd.read_excel("20660-23_final_processed.xlsx")

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
with open('DNPM-BN-0001.json', 'w') as file:
    json.dump(
        DNPM_out, 
        file, 
        indent=4, 
        sort_keys=True
    )

# End    
###############################################################################
