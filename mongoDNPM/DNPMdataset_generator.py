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
file_dnpm_data_template="PATIENT_ID1_DNPM_DatenmodellV2.0.xlsx"

# Get dnpm_data
with open (file_dnpm_data_template, "rb") as f:
    dnpm_data = pd.read_excel(f,\
                         sheet_name=None, engine = "openpyxl")
        
        
# get "Einfache Varianten"
processed_variant_data = pd.read_excel("E22533-23_final_processed.xlsx")

# get relevant data for report

# get idx of reported
idx_of_variants_to_report = processed_variant_data.index\
                           [processed_variant_data["Submit"] == "x"].tolist()

variants_to_report = processed_variant_data.loc[idx_of_variants_to_report, :]

# get all necessary columns
final_data = variants_to_report[["Chromosome",
                                 "SYMBOL",
                                 "NM_v",
                                 "Exon Number",
                                "Position",
                                "End Position",
                                "Reference",
                                "Allele",
                                "HGVSc_x",
                                "HGVS_PROTEIN",
                                "Coverage",
                                "Frequency",
                                "name dbsnp_v151_ensembl_hg38_no_alt_analysis_set",
                                "CLIN_SIG"]]

# rename to DNPMdataset names
final_data = final_data.rename(columns=\
           {"Chromosome": "Chromosom",
            "SYMBOL": "Gen",
            "NM_v": "Transcript",
            "Exon Number": "Exon",
            "Reference": "Ref",
            "Allele": "Alt",
            "HGVSc_x": "cDNA Change",
            "HGVS_PROTEIN": "Amino Acid Change",
            "Coverage": "Read Depth",
            "Frequency": "Allelic Frequency", 
            "name dbsnp_v151_ensembl_hg38_no_alt_analysis_set": "dbSNP ID",
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
DNPM_out = dnpm_data["NGS-Bericht"].to_dict(orient='records')

# convert dataframes into dictionaries
dnpm_data_dict = {
    key: dnpm_data[key].to_dict(orient='records') 
    for key in dnpm_data.keys()
}

# add variant information
DNPM_out[0]["Einfache Varianten"] = dnpm_data_dict["Einfache Variante"]

def unlist_func(DNPM_out):
    
    for i in DNPM_out:
        return i
    
DNPM_out = unlist_func(DNPM_out)

# write to disk
with open('data_dict.json', 'w') as file:
    json.dump(
        DNPM_out, 
        file, 
        indent=4, 
        sort_keys=True
    )

# End    
###############################################################################
