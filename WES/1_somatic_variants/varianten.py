# -*- coding: utf-8 -*-
"""
Created on Sat Dec  9 12:16:51 2023

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
parser.add_argument("-q", "--vep_file", type=str)
parser.add_argument("-c", "--clc_file", type=str)
#parser.add_argument("-l", "--gene_list", type=str)
parser.add_argument("-t", "--transcript_list", type=str)  
parser.add_argument("-o", "--output", type=str)   
parser.add_argument("-tmb", "--tmb_output", type=str)              
args = parser.parse_args()

#Fall = args.Fall
#vep_file = "X:/PAT-Sequenzer/WES-Pilot/WES_Pilot_input_all/snv_vep/14_T_e.txt"
#clc_file = "X:/PAT-Sequenzer/WES-Pilot/WES_Pilot_final/somatic_variants/txt/Final_filtered_variants-Fall_14_T_1 (paired), Fall_14_N_1 (paired).txt"
#Genliste_somatisch = args.Genliste_somatisch
#Genliste_mitTr_somatisch = args.Genliste_mitTr_somatisch
#output = "./"

# hg 38
# 1.0 VEP data
VEP_data = pd.read_csv(args.vep_file, delimiter="\t")  
#VEP_data = pd.read_csv(vep_file, delimiter="\t")
# 2.0 CLC variant track data
CLC_variant_track_data = pd.read_csv(args.clc_file, delimiter="\t") 
#CLC_variant_track_data = pd.read_csv(clc_file, delimiter="\t") 

#------------------------------------------------------------------------------

# 2.1 rename Region to position
#VEP_data = VEP_data.rename\
#                        (columns={"POS": "Position"}) # not possible here!
                        
CLC_variant_track_data = CLC_variant_track_data.rename\
                        (columns={"Region": "Position"})
                        
# 2.2 add new column End Position
CLC_variant_track_data.insert(loc=2, column="End Position", value="")

# 2.3 rename CHROM to Chromosome
#VEP_data = VEP_data.rename\
#                        (columns={"CHROM": "Chromosome"}) # not possible here!

# CLC data
for i in range(len(CLC_variant_track_data["Position"])):
    
    if ".." not in CLC_variant_track_data["Position"][i] and \
        "^" not in CLC_variant_track_data["Position"][i]:
            tmp_value_snv = CLC_variant_track_data["Position"][i]
            CLC_variant_track_data.loc\
            [CLC_variant_track_data.index[i], "End Position"] = tmp_value_snv
            
    if ".." in CLC_variant_track_data["Position"][i]:
            tmp_value_del = CLC_variant_track_data["Position"][i]
            tmp_value_del_start = tmp_value_del.split("..")[0]
            tmp_value_del_end = tmp_value_del.split("..")[1]
            CLC_variant_track_data.loc\
            [CLC_variant_track_data.index[i], "Position"] =  tmp_value_del_start
            CLC_variant_track_data.loc\
            [CLC_variant_track_data.index[i], "End Position"] = tmp_value_del_end
            
    if "^" in CLC_variant_track_data["Position"][i]:
            tmp_value_ins = CLC_variant_track_data["Position"][i]
            tmp_value_ins_start = tmp_value_ins.split("^")[0]
            tmp_value_ins_end = tmp_value_ins.split("^")[1]
            CLC_variant_track_data.loc\
            [CLC_variant_track_data.index[i], "Position"] = tmp_value_ins_start
            CLC_variant_track_data.loc\
            [CLC_variant_track_data.index[i], "End Position"] = tmp_value_ins_end

#------------------------------------------------------------------------------

# Filtering of valid variants

# rename column Reference allele to ReferenceAllele
CLC_variant_track_data = CLC_variant_track_data.rename(columns={"Reference allele": "ReferenceAllele"})
group_reference_allele = CLC_variant_track_data.groupby(CLC_variant_track_data.ReferenceAllele)
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
    clc_data_ReferenceAllele_NO[columns_to_filter[i]] = variants_lst_sheets_temp4
 
# Filter data with QUAL >= 150
clc_data_filtered = clc_data_ReferenceAllele_NO\
[clc_data_ReferenceAllele_NO["QUAL"] >= 150]

# discarded 1
#clc_data_discarded_1 = clc_data_ReferenceAllele_NO\
#[clc_data_ReferenceAllele_NO["QUAL"] < 150]
 
# Filter data with Av quality >= 25
clc_data_filtered = clc_data_filtered\
[clc_data_filtered["Average quality"] >= 25]

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
 
clc_data_filtered = clc_data_filtered.reset_index()

#------------------------------------------------------------------------------

# start here TMB
get_singe_mut = clc_data_filtered["Position"].drop_duplicates()
get_singe_mut_idx = get_singe_mut.index.tolist()

TMB_df = clc_data_filtered.loc[get_singe_mut_idx, :]

# get snv
TMB_snv = TMB_df[TMB_df["Type"] == "SNV"]
TMB_snv = TMB_snv[TMB_snv["Non-synonymous"] != "No"]
TMB_snv = TMB_snv[TMB_snv["Non-synonymous"] == "Yes"]
TMB_snv = TMB_snv[TMB_snv["Coverage"] >= 100]
#TMB_snv = TMB_snv[TMB_snv["Frequency"] >= 5]

TMB_del = TMB_df[TMB_df["Type"] == "Deletion"]
TMB_del = TMB_del [TMB_del ["Non-synonymous"] != "No"]
TMB_del  = TMB_del [TMB_del ["Non-synonymous"] == "Yes"]
TMB_del = TMB_del[TMB_del["Coverage"] >= 100]
TMB_in = TMB_df[TMB_df["Type"] == "Insertion"]
TMB_in = TMB_in[TMB_in["Non-synonymous"] != "No"]
TMB_in = TMB_in[TMB_in["Non-synonymous"] == "Yes"]
TMB_in = TMB_in[TMB_in["Coverage"] >= 100]

TMB_snv_final = len(TMB_snv)
TMB_indel = len(TMB_snv) + len(TMB_del) + len(TMB_in)

#Regionsgroesse_MB = 33.94
Regionsgroesse_MB = 30.16

# make dataframe
TMB = pd.DataFrame()

header_col = ["Anzahl TMB Mut. missense", "Anzahl TMB Mut. Missense + InDel", "Regionsgroesse [Mb]",\
              "TMB Missense", \
              "TMB Missense + InDel"]

for col in header_col:
    TMB [col] = ""
    
TMB.loc[0, "Anzahl TMB Mut. missense"] = TMB_snv_final
TMB.loc[0, "Anzahl TMB Mut. Missense + InDel"] = TMB_indel
TMB.loc[0, "Regionsgroesse [Mb]"] = Regionsgroesse_MB
TMB.loc[0, "TMB Missense"] = round(TMB_snv_final/Regionsgroesse_MB, 2)
TMB.loc[0, "TMB Missense + InDel"] = round(TMB_indel/Regionsgroesse_MB, 2)
    
TMB.to_csv(args.tmb_output, index=False) 

#------------------------------------------------------------------------------

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
#transcript_list = "X:/PAT-Sequenzer/WES-Pilot/WES_Pilot_input_all/Genliste_mitTr_somatisch.xlsx"
transcript_list = args.transcript_list
RefSeq_NM = pd.read_excel(transcript_list)
RefSeq_NM_lst = RefSeq_NM["NM_RefSeq_final"].values.tolist()
# ok bis hier1
#------------------------------------------------------------------------------

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
    
   
# check for duplicate values in column
#duplicate_values = pre_final_data["NM"].duplicated()
#print(duplicate_values)

# check for duplicate values in  column
#duplicate_values = clc_data_filtered_dropXM["Position"].duplicated()
#print(duplicate_values)

# remove duplicate values in ProductID column
#df = pre_final_data.drop_duplicates(subset=['NM'], keep='first')
#------------------------------------------------------------------------------

# 4.0 get gene list
#gene_list = "X:/PAT-Sequenzer/WES-Pilot/WES_Pilot_input_all/Genliste_somatisch.csv"
#gene_list_somatic = pd.read_csv(gene_list, header=None)
#gene_list_somatic = gene_list_somatic.rename(columns={0: "Gene"})
#gene_list_somatic = gene_list_somatic["Gene"].values.tolist()

# 5.0 filter after gene lst somatic
#clc_data_filtered_genelist = clc_data_filtered
#true_idx = []
#false_idx = []

#for gene in range(len(clc_data_filtered_genelist["Homo_sapiens_refseq_GRCh38.p13_no_alt_analysis_set_Genes"])):
#    if clc_data_filtered_genelist["Homo_sapiens_refseq_GRCh38.p13_no_alt_analysis_set_Genes"][gene] in gene_list_somatic:
#        true_idx.append(gene)
#    if clc_data_filtered_genelist["Homo_sapiens_refseq_GRCh38.p13_no_alt_analysis_set_Genes"][gene] not in gene_list_somatic:
#        false_idx.append(gene)

# 6.0 pre_results
#pre_final = clc_data_filtered_genelist.loc[true_idx, :]
#pre_discarded = clc_data_filtered_genelist.loc[false_idx, :]

#pre_final = merged.loc[true_idx, :]
#pre_discarded = merged.loc[false_idx, :]

#------------------------------------------------------------------------------

# prcessing vep data
                        
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
    
#------------------------------------------------------------------------------

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


merged = pd.merge(pre_final_data, VEP_data[["Chromosome", "Position",\
                 "NM-Nummer", "HGVSc", "HGVSp", "SYMBOL"]], \
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


        
final_variants = merged.sort_values(by=["Frequency"], ascending=False)                

final_variants.to_excel(args.output, \
                             index = False, \
                             engine= None) 
  
