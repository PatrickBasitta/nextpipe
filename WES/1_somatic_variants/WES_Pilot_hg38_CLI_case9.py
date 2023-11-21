# -*- coding: utf-8 -*-
"""
Created on Thu Oct  5 09:34:34 2023

Pyhton script for merging/joining csv data obtained from CLC QCI and CLC variant
track, to filter accodingly and convert results to csv WES-Pilot format

ONLY USE FOR WES-PILOT

@author: Patrick Basitta
"""
import pandas as pd
import numpy as np
import csv
import os
import re
import argparse

#---------------------------------------
# Using argparse for positinal arguments
#---------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-n", "--number", type=str)
parser.add_argument("-q", "--qci_file", type=str)
parser.add_argument("-c", "--clc_file", type=str)
parser.add_argument("-l", "--gene_list", type=str)
parser.add_argument("-t", "--transcript_list", type=str)  
parser.add_argument("-o", "--output", type=str)                
args = parser.parse_args()

#Fall = args.Fall
#QCI_filename = args.QCI_file
#CLC_filename = args.CLC_file
#Genliste_somatisch = args.Genliste_somatisch
#Genliste_mitTr_somatisch = args.Genliste_mitTr_somatisch
#Fall_output = args.Fall_output
# hg 38
# 1.0 QCI data
QCI_data = pd.read_csv(args.qci_file, delimiter="\t")  

# 2.0 CLC variant track data
CLC_variant_track_data = pd.read_csv(args.clc_file, delimiter="\t") 

#------------------------------------------------------------------------------

# 2.1 rename Region to position
QCI_data = QCI_data.rename\
                        (columns={"Region": "Position"})

CLC_variant_track_data = CLC_variant_track_data.rename\
                        (columns={"Region": "Position"})

# 2.2 add new column End Position 
QCI_data.insert(loc=2, column="End Position", value="")
CLC_variant_track_data.insert(loc=2, column="End Position", value="")

# 2.3 add values SNV, Indels to columns Position and End Position 
# QCI data
for i in range(len(QCI_data["Position"])):
    
    if ".." not in QCI_data["Position"][i] and \
        "^" not in QCI_data["Position"][i]:
            tmp_value_snv = QCI_data["Position"][i]
            QCI_data.loc\
            [QCI_data.index[i], "End Position"] = tmp_value_snv
            
    if ".." in QCI_data["Position"][i]:
            tmp_value_del = QCI_data["Position"][i]
            tmp_value_del_start = tmp_value_del.split("..")[0]
            tmp_value_del_end = tmp_value_del.split("..")[1]
            QCI_data.loc\
            [QCI_data.index[i], "Position"] =  tmp_value_del_start
            QCI_data.loc\
            [QCI_data.index[i], "End Position"] = tmp_value_del_end
            
    if "^" in QCI_data["Position"][i]:
            tmp_value_ins = QCI_data["Position"][i]
            tmp_value_ins_start = tmp_value_ins.split("^")[0]
            tmp_value_ins_end = tmp_value_ins.split("^")[1]
            QCI_data.loc\
            [QCI_data.index[i], "Position"] = tmp_value_ins_start
            QCI_data.loc\
            [QCI_data.index[i], "End Position"] = tmp_value_ins_end

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
       
# 3.0 merge QCI and CLC dfs
QCI_data["Chromosome"] = QCI_data["Chromosome"].astype(str)
CLC_variant_track_data["Chromosome"] = CLC_variant_track_data["Chromosome"].astype(str)

QCI_data["Position"] = QCI_data["Position"].astype(int)
CLC_variant_track_data["Position"] = CLC_variant_track_data["Position"].astype(int)

QCI_data["End Position"] = QCI_data["End Position"].astype(int)
CLC_variant_track_data["End Position"] = CLC_variant_track_data["End Position"].astype(int)

merged = pd.merge(QCI_data, CLC_variant_track_data[["Chromosome", "Position",\
                 "End Position", "Coverage", "Frequency", \
                 "Forward/reverse balance", \
                 "BaseQRankSum", "Read position test probability", \
                 "Read direction test probability", "Homopolymer length", \
                 "Homopolymer", "QUAL", "Non-synonymous"]], \
                 on = ["Chromosome", "Position", "End Position"], how = "left")
  
# 4.0 get gene list
gene_list_somatic = pd.read_csv(args.gene_list, header=None)
gene_list_somatic = gene_list_somatic.rename(columns={0: "Gene"})
gene_list_somatic = gene_list_somatic["Gene"].values.tolist()

# 5.0 filter after gene lst somatic
true_idx = []
false_idx = []

for gene in range(len(merged["GENE_SYMBOL"])):
    if merged["GENE_SYMBOL"][gene] in gene_list_somatic:
        true_idx.append(gene)
    if merged["GENE_SYMBOL"][gene] not in gene_list_somatic:
        false_idx.append(gene)
        
# 6.0 pre_results
pre_final = merged.loc[true_idx, :]
pre_discarded = merged.loc[false_idx, :]

#pre_final["GENE_SYMBOL"]

# 6.1 filter RefSeq
RefSeq_NM = pd.read_excel(args.transcript_list)
RefSeq_NM_lst = RefSeq_NM["NM_RefSeq_final"].values.tolist()

# Reset indices
pre_final = pre_final.reset_index()

NM_idx = []
for i in range(len(pre_final["TRANSCRIPT_ID"])):
    if pre_final["TRANSCRIPT_ID"][i].split(".")[0] in RefSeq_NM_lst:
        NM_idx.append(i)
        
pre_final_data = pre_final.loc[NM_idx, :]


# 6.2 get Pathogenic, Likely_Pathogenic and risk factor
#CLI_idx = []
#ING_idx = []
#CLNSIG_idx = []
#missed_gene_idx = []
#classification = ["Pathogenic", "Likely_Pathogenic", "Likely_pathogenic", "risk_factor"]
#for i in range(len(pre_final["CLI_ASSESSMENT"])):
#    if pre_final["CLI_ASSESSMENT"][i] in classification:
#        CLI_idx.append(i)
#    if pre_final["ING_CLASSIFICATION"][i] in classification:
#        ING_idx.append(i)
#    if pre_final["CLNSIG clinvar_20210828_hg38_no_alt_analysis_set"][i] in classification:
#        CLNSIG_idx.append(i)
 
#combine_lst = set(CLI_idx + ING_idx + CLNSIG_idx)  

#combine_lst = list(combine_lst)

#PLP_pre_final = pre_final.loc[combine_lst, :]

#gene_lst = pre_final_data["GENE_SYMBOL"].tolist()
#NM_used = pre_final_data["TRANSCRIPT_ID"].tolist()

#for gene in range(len(PLP_pre_final["GENE_SYMBOL"])):
#    if PLP_pre_final["GENE_SYMBOL"].iloc[gene] not in gene_lst and \
#        PLP_pre_final["TRANSCRIPT_ID"].iloc[gene] not in NM_used: 
#            missed_gene_idx.append(gene)
            
# concate missed gen with pre_final_data df

#PLP_pre_final_missed = PLP_pre_final.iloc[missed_gene_idx]

#pre_final_data = pd.concat([pre_final_data, PLP_pre_final_missed])

# 7.0 filter after QUAL and for/rev balance  

# retyping
#type(pre_final["QUAL"].iloc[0])    # str 
pre_final_temp1_Q = pre_final_data["QUAL_y"].apply(lambda x: x.strip("'"))
pre_final_temp2_Q = pre_final_temp1_Q.apply(lambda x: x.replace(",","."))
pre_final_temp3_Q = pre_final_temp2_Q.to_frame()
pre_final_data["QUAL_y"] = pre_final_temp3_Q.astype("float")

# step 7.1 filter QUAL                            
#tmp_final_sorted_Q = pre_final[pre_final["QUAL"] >= 200] 
tmp_pre_final_Q = pre_final_data\
    [pre_final_data["QUAL_y"]\
     >= 200]

# discarded variatns #1
#tmp_discarded_Q = pre_final[pre_final["QUAL"] < 200]
tmp_discarded_Q = pre_final_data\
    [pre_final_data["QUAL_y"]\
     < 200]

# step 7.2 for/rev balance (not available after remap)
#tmp_pre_final_Q_rfb = tmp_pre_final_Q\
#    [tmp_pre_final_Q["Forward/reverse balance"] != "0"]

#tmp_discarded_rfb = tmp_discarded_Q\
#    [tmp_discarded_Q["Forward/reverse balance"] == "0"]
    
# 8.0 sort frequency (not necessary)
#pre_final_temp1_QF = tmp_pre_final_Q["ING_AF"].apply(lambda x: x.strip("'"))
#pre_final_temp2_QF = pre_final_temp1_QF.apply(lambda x: x.replace(",","."))
#pre_final_temp3_QF= pre_final_temp2_QF.to_frame()
#tmp_pre_final_Q["ING_AF"] = pre_final_temp3_QF.astype("float")

# 9.0 
final_variants = tmp_pre_final_Q 

#final_variants_col = final_variants[["Chromosome",
#                                "Position",
#                                "End Position",
#                                "Reference",
#                                "Allele",
#                                "Coverage_x",
#                                "ING_AF",
#                                "BaseQRankSum",
#                                "Read position test probability",
#                                "Read direction test probability", 
#                                "Homopolymer length", 
#                                "Homopolymer", 
#                                "QUAL_y",
#                                "FILTER",
#                                "TRANSLATION_IMPACT",
#                                "GENE_SYMBOL",
#                                "TRANSCRIPT_ID",
#                                "HGVS_PROTEIN",
#                                "HGVS_TRANSCRIPT",
#                                "Non-synonymous",
#                                "GENE_REGION",
#                                "ING_IA",
#                                "CLI_ASSESSMENT",
#                                "ING_CLASSIFICATION",
#                                "CLNSIG clinvar_20210828_hg38_no_alt_analysis_set",
#                                "SIFT_FUNCTION",
#                                "POLYPHEN_FUNCTION",
#                                "1000G_AF"]]

# 10.0 write to xlsx
final_variants.to_excel(args.output, \
                             index = False, \
                             engine= None)


