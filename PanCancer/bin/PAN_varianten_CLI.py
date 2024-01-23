#!/usr/bin/env python
"""
Created on Mon Jan 15 14:48:10 2024

@author: PatrickBasitta
"""

import pandas as pd
import numpy as np
import csv
import os
import re
import argparse
import json
from itertools import zip_longest

#---------------------------------------
# Using argparse for positinal arguments
#---------------------------------------
parser = argparse.ArgumentParser()
#parser.add_argument("-n", "--number", type=str)
parser.add_argument("-q", "--vep_PAN_file", type=str)
parser.add_argument("-c", "--clc_PAN_file", type=str)
#parser.add_argument("-l", "--gene_list", type=str)
parser.add_argument("-t", "--transcript_PAN_list", type=str)  
parser.add_argument("-o", "--outfile", type=str)   
#parser.add_argument("-tmb", "--tmb_output", type=str)              
args = parser.parse_args()

#--------------------------------------
# Getting CLC PAN DATA and format
#-------------------------------------
#clc_PAN_file = "/data2/basitta/pan_data/*.txt"
clc_PAN_file = args.clc_PAN_file
CLC_variant_track_data_PAN = pd.read_csv(clc_PAN_file, delimiter=";", encoding="ISO-8859-1") 
#CLC_variant_track_data_PAN = pd.read_csv(clc_PAN_file)
# Rename Region to position
CLC_variant_track_data_PAN = CLC_variant_track_data_PAN.rename\
                        (columns={"Region": "Position"})
                        
# Add new column End Position
CLC_variant_track_data_PAN.insert(loc=2, column="End Position", value="")

# CLC PAN data
for i in range(len(CLC_variant_track_data_PAN["Position"])):
    
    if ".." not in CLC_variant_track_data_PAN["Position"][i] and \
        "^" not in CLC_variant_track_data_PAN["Position"][i]:
            tmp_value_snv = CLC_variant_track_data_PAN["Position"][i]
            CLC_variant_track_data_PAN.loc\
            [CLC_variant_track_data_PAN.index[i], "End Position"] = tmp_value_snv
            
    if ".." in CLC_variant_track_data_PAN["Position"][i]:
            tmp_value_del = CLC_variant_track_data_PAN["Position"][i]
            tmp_value_del_start = tmp_value_del.split("..")[0]
            tmp_value_del_end = tmp_value_del.split("..")[1]
            CLC_variant_track_data_PAN.loc\
            [CLC_variant_track_data_PAN.index[i], "Position"] =  tmp_value_del_start
            CLC_variant_track_data_PAN.loc\
            [CLC_variant_track_data_PAN.index[i], "End Position"] = tmp_value_del_end
            
    if "^" in CLC_variant_track_data_PAN["Position"][i]:
            tmp_value_ins = CLC_variant_track_data_PAN["Position"][i]
            tmp_value_ins_start = tmp_value_ins.split("^")[0]
            tmp_value_ins_end = tmp_value_ins.split("^")[1]
            CLC_variant_track_data_PAN.loc\
            [CLC_variant_track_data_PAN.index[i], "Position"] = tmp_value_ins_start
            CLC_variant_track_data_PAN.loc\
            [CLC_variant_track_data_PAN.index[i], "End Position"] = tmp_value_ins_end
            
# Filtering of valid variants

# rename column Reference allele to ReferenceAllele
CLC_variant_track_data_PAN = CLC_variant_track_data_PAN.rename(columns={"Reference allele": "ReferenceAllele"})
group_reference_allele = CLC_variant_track_data_PAN.groupby(CLC_variant_track_data_PAN.ReferenceAllele)
clc_data_ReferenceAllele_NO =  group_reference_allele.get_group("No")
 
columns_to_filter = ["Frequency", "QUAL", "Forward/reverse balance", \
                     "Average quality",\
                     "Read position test probability", \
                     "Read direction test probability"]
     
# Convert string number in float of columns to filter   
for i in range(len(columns_to_filter)):
 
    variants_lst_sheets_temp1 = clc_data_ReferenceAllele_NO\
        [columns_to_filter[i]].apply(lambda x: x.strip("'"))
    variants_lst_sheets_temp2 = variants_lst_sheets_temp1.apply\
                                             (lambda x: x.replace(",","."))
    variants_lst_sheets_temp3 = variants_lst_sheets_temp2.to_frame()
    variants_lst_sheets_temp4 = variants_lst_sheets_temp3.astype("float")
    clc_data_ReferenceAllele_NO.loc[:,columns_to_filter[i]] = variants_lst_sheets_temp4
    
#type(clc_data_ReferenceAllele_NO.loc[1,"Frequency"])
    
# Filter data with QUAL >= 150
clc_data_filtered = clc_data_ReferenceAllele_NO\
[clc_data_ReferenceAllele_NO["QUAL"] >= 150]

# discarded 1
#clc_data_discarded_1 = clc_data_ReferenceAllele_NO\
#[clc_data_ReferenceAllele_NO["QUAL"] < 150]
 
# Filter data with Av quality >= 25
clc_data_filtered = clc_data_filtered\
[clc_data_filtered["Average quality"] >= 35]

# discarded 2
#clc_data_discarded_2 = clc_data_filtered\
#[clc_data_filtered["Average quality"] < 25]
 
# Filter data with count >= 2
clc_data_filtered = clc_data_filtered\
[clc_data_filtered["Count"] >= 2]

# discarded 3
#clc_data_discarded_3 = clc_data_filtered\
#[clc_data_filtered["Count"] < 2]
 
# Filter data with Frequency >= 5
clc_data_filtered = clc_data_filtered\
[clc_data_filtered["Frequency"] >= 5]

# discarded 4
#clc_data_discarded_4 = clc_data_filtered\
#[clc_data_filtered["Frequency"] < 5]
 
# Filter data with Forward/reverse balance > 0
clc_data_filtered = clc_data_filtered\
[clc_data_filtered["Forward/reverse balance"] > 0]
 
# Filter data with Read position test probability >= 0.000001
clc_data_filtered = clc_data_filtered\
[clc_data_filtered["Read position test probability"] >= 0.000001]
 
# Filter data with Read direction test probability >= 0.000001
clc_data_filtered = clc_data_filtered\
[clc_data_filtered["Read direction test probability"] >= 0.000001]

# filter non-synonymous???
 
clc_data_filtered = clc_data_filtered.reset_index()

# "explode" coding region change and set NM and HGVSc column           
clc_data_filtered = clc_data_filtered.rename(columns=\
           {"Coding region change": "Coding_region_change"})

clc_data_filtered = (
    clc_data_filtered.assign(Coding_region_change=clc_data_filtered\
                                  ["Coding_region_change"].str.split("; "))
        .explode("Coding_region_change")
        .reset_index(drop=True)
)
    
for index, NM_tr in enumerate(clc_data_filtered["Coding_region_change"]):
    
    if pd.isna(NM_tr):
        clc_data_filtered.loc[i, "NM"] = clc_data_filtered.loc[i, "NM"]
        clc_data_filtered.loc[index, "HGVSc"] = clc_data_filtered.loc[index, "HGVSc"]
        
    else:
        clc_data_filtered.loc[index, "NM"] = NM_tr.split(":")[0]
        clc_data_filtered.loc[index, "HGVSc"] = NM_tr.split(":")[1]
        
#------------------------------------------------------------------------------

# Drop XM transcripts and filter after transcriptlist !not necessary
#clc_data_filtered_dropXM = clc_data_filtered[clc_data_filtered\
#                           ["NM"].str.contains("XM*") == False]
clc_data_filtered_dropNa = clc_data_filtered[clc_data_filtered\
                          ["NM"].str.contains("nan") == False]
# clc_data_filtered_XM = clc_data_filtered[clc_data_filtered\
#                            ["NM"].str.contains("XM*") == True]
#------------------------------------------------------------------------------

# filter RefSeq
clc_data_filtered_dropXM = clc_data_filtered_dropNa
#transcript_list = "X:/PAT-Sequenzer/PanCancer_test/finale_aktuelle_pancancer_transcript_lst.xlsx"
transcript_list = args.transcript_PAN_list
#transcript_list = args.transcript_list
RefSeq_NM = pd.read_excel(transcript_list)
RefSeq_NM_lst = RefSeq_NM["NM_RefSeq_final"].values.tolist()
# ok bis hier1
# in list with version number or in clc_data_filtered_dropXM!!!

# Reset indices
clc_data_filtered_dropXM = clc_data_filtered_dropXM.reset_index()
clc_data_filtered_dropXM["NM"] = clc_data_filtered_dropXM["NM"].astype(str)
NM_idx = []
for i in range(len(clc_data_filtered_dropXM["NM"])):
    if clc_data_filtered_dropXM["NM"][i].split(".")[0] in RefSeq_NM_lst:
        NM_idx.append(i)

pre_final_data = clc_data_filtered_dropXM.loc[NM_idx, :]
# ok bis hier2

index_variants = pre_final_data.index.tolist()
for i in range(len(pre_final_data["NM"])):
    tmp_NM_name = pre_final_data["NM"][index_variants[i]].split(".")
    pre_final_data.loc[[index_variants[i]], "NM"] =  tmp_NM_name[0]

# before vep annotation with vep
# final_variants = pre_final_data.sort_values(by=["Frequency"], ascending=False)                

# final_variants.to_excel("test_output.xlsx", \
#                             index = False, \
#                             engine= None) 

# prcessing vep data
#vep_file = "X:/PAT-Sequenzer/PanCancer_test/*.txt"
vep_file = args.vep_PAN_file
VEP_data = pd.read_csv(vep_file, delimiter="\t")                        
# add new column Chromosome and Position
VEP_data.insert(loc=1, column="Chromosome", value="")
VEP_data.insert(loc=1, column="Position", value="")

# VEP data
for i in range(len(VEP_data["Location"])):
    tmp_split = VEP_data["Location"][i].split("-")[0]
    tmp_split = tmp_split.split(":")
    VEP_data.loc[i, "Chromosome"] =  tmp_split[0]
    VEP_data.loc[i, "Position"] =  tmp_split[1]

for i in range(len(VEP_data["Feature"])):
    tmp_NM = VEP_data["Feature"][i].split(".")
    VEP_data.loc[i, "Feature"] =  tmp_NM[0]

# 3.0 merge VEP and CLC dfs
pre_final_data["Chromosome"] = pre_final_data["Chromosome"]
VEP_data["Chromosome"] = VEP_data["Chromosome"]

VEP_data["Chromosome"] = VEP_data["Chromosome"].astype(str)
pre_final_data["Chromosome"] = pre_final_data["Chromosome"].astype(str)

VEP_data["Position"] = VEP_data["Position"].astype(int)
pre_final_data["Position"] = pre_final_data["Position"].astype(int)
                 
pre_final_data = pre_final_data.rename\
                        (columns={"NM": "NM-Nummer"})

VEP_data = VEP_data.rename\
                        (columns={"Feature": "NM-Nummer"})
                        
VEP_data["NM-Nummer"] = VEP_data["NM-Nummer"].astype(str)
pre_final_data["NM-Nummer"] = pre_final_data["NM-Nummer"].astype(str)

# add here more columns from vep
merged = pd.merge(pre_final_data, VEP_data[["Chromosome", "Position",\
                 "NM-Nummer", "HGVSc", "HGVSp", "SYMBOL", "AF", "MAX_AF", \
                 "gnomADe_AF", "gnomADg_AF", "SIFT", "PolyPhen", "PUBMED" ]], \
                 on = ["Chromosome", "Position", "NM-Nummer"], how = "left")
    
merged = merged.rename(columns={"NM-Nummer": "TRANSCRIPT_ID"})

index_variants_merged = merged.index.tolist()
for i in range(len(merged["HGVSp"])):
    if pd.isna(merged["HGVSp"][index_variants_merged[i]]):
        merged.loc[[index_variants_merged[i]], "HGVS_PROTEIN"] = merged["HGVSp"][index_variants_merged[i]]
    else:
        tmp_p_name = merged["HGVSp"][index_variants_merged[i]].split(":")
        if len(tmp_p_name) >= 2:
            merged.loc[[index_variants_merged[i]], "HGVS_PROTEIN"] =  tmp_p_name[1]

# get comprehensive output format

processed_data_final = merged[["Chromosome", "Position", "End Position", \
                       "Reference", "Allele", "Count", "Coverage", \
                        "Frequency", "QUAL", "Forward/reverse balance", \
                        "Average quality","Read position test probability", \
                        "Read direction test probability", "BaseQRankSum",
                        "Homopolymer", "Homopolymer length", \
                        "Count (singleton UMI)", "Count (big UMI)", \
                        "Proportion (singleton UMIs)", "SYMBOL", \
                        "TRANSCRIPT_ID", "HGVSc_x", "HGVS_PROTEIN", "Exon Number", \
                        "func dbsnp_v151_ensembl_hg38_no_alt_analysis_set", \
                        "name dbsnp_v151_ensembl_hg38_no_alt_analysis_set", \
                        "CLNSIG clinvar_20220730_hg38_no_alt_analysis_set", \
                        "CLNREVSTAT clinvar_20220730_hg38_no_alt_analysis_set", \
                        "AF_EXAC clinvar_20220730_hg38_no_alt_analysis_set", \
                        "AF_TGP clinvar_20220730_hg38_no_alt_analysis_set", \
                        "AF", "MAX_AF", "gnomADe_AF", "gnomADg_AF", "SIFT", \
                        "PolyPhen", "PUBMED"]]

 # round AF to max 2 decimals (since AF number are str in excel, one has to set these values to int before run this script!
processed_data_final.loc[:,"Frequency"]  = processed_data_final\
                                   ["Frequency"].apply(lambda x: round(x,2))

     
final_variants = processed_data_final.sort_values(by=["Frequency"], ascending=False)                

#final_variants.to_excel("final_output_with_aktuelle_PanCancerliste_AF.xlsx", \
#                             index = False, \
#                             engine= None) 

final_variants.to_excel(args.outfile, \
                            index = False, \
                             engine= None) 






































