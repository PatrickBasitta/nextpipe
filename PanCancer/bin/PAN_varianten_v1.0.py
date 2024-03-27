#!/usr/bin/env python
"""
Created on Mon Jan 15 14:48:10 2024

Script to process variant output obtained from CLC PanCancer Workflow and
ENSEMBL VEP TOOL

Input: clc_PAN_file as .csv
       vep_PAN_file as .txt
       transcript_PAN_list as .xlsx
       variantDBi as .xlsx

output: user specified name; file extension: .xlsx

@author: PatrickBasitta
"""

import pandas as pd
import argparse
import functions_pan.process_variantlist_PAN_utils as fp
from datetime import datetime

#---------------------------------------
# Using argparse for positinal arguments
#---------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-c", "--clc", type=str)
parser.add_argument("-v", "--vep", type=str)
parser.add_argument("-t", "--transcripts", type=str)
parser.add_argument("-D", "--variant_DBi", type=str)
parser.add_argument("-o", "--outfile", type=str)
parser.add_argument("-rv", "--removed_variants", type=str)
parser.add_argument("-rd", "--removed_duplicates", type=str)
args = parser.parse_args()

#--------------------------------------------------
# Getting current date and time for log information
#--------------------------------------------------
date_time_now = datetime.now()
# dd/mm/YY H:M:S
dt_string = date_time_now.strftime("%d/%m/%Y %H:%M:%S")
print("Start:", dt_string)
# Processed file
print("Input_CLC_file:", args.clc)
print("Input_VEP_file:", args.vep)
print("Output_file_1:", args.outfile)
print("Output_file_2:", args.removed_variants)
print("Output_file_3:", args.removed_duplicates)
# Script used
print("Script: PAN_varianten_v1.0.py")

#-----------------
# Get CLC_PAN_data
#-----------------
# clc_PAN_file = ".csv"
clc_PAN_file = args.clc
CLC_variant_track_data_PAN = pd.read_csv(clc_PAN_file, delimiter=";",\
                                       encoding="ISO-8859-1")
# CLC_variant_track_data_PAN = pd.read_csv(clc_PAN_file, delimiter=",")

#--------------------------------------
# CLC_PAN_data - adjust region_position
#--------------------------------------
CLC_variant_track_data_PAN = fp.adjust_region_position(CLC_variant_track_data_PAN)

#------------------------------------------------------          
# Filtering of variants not equal to "Reference allele"
# Rename column "Reference allele" to "ReferenceAllele"
#------------------------------------------------------
CLC_variant_track_data_PAN = CLC_variant_track_data_PAN.rename\
                              (columns={"Reference allele": "ReferenceAllele"})
group_reference_allele = CLC_variant_track_data_PAN.groupby\
                                   (CLC_variant_track_data_PAN.ReferenceAllele)
clc_data_ReferenceAllele_NO =  group_reference_allele.get_group("No")

#----------------------
# Column list to filter
#----------------------
columns_to_filter = ["Frequency", "QUAL", "Forward/reverse balance", \
                     "Average quality",\
                     "Read position test probability", \
                     "Read direction test probability"]

#----------------------------------------------------   
# Convert string number in float of columns_to_filter
#----------------------------------------------------   
clc_data_ReferenceAllele_NO = fp.convert_coltype_str_to_float\
    (columns_to_filter, clc_data_ReferenceAllele_NO)
       
#-----------------------------
# Filter data with QUAL >= 150
#-----------------------------
clc_data_filtered = clc_data_ReferenceAllele_NO\
[clc_data_ReferenceAllele_NO["QUAL"] >= 150]

# discarded variants 1 (keep discarded data)
clc_data_discarded_1 = clc_data_ReferenceAllele_NO\
[clc_data_ReferenceAllele_NO["QUAL"] < 150]

#---------------------------------- 
# Filter data with Av quality >= 35
#----------------------------------
clc_data_filtered_aq = clc_data_filtered\
[clc_data_filtered["Average quality"] >= 35]

# discarded variants 2 (keep discared data)
clc_data_discarded_2 = clc_data_filtered\
[clc_data_filtered["Average quality"] < 25]

#----------------------------
# Filter data with count >= 2
#----------------------------
clc_data_filtered_c = clc_data_filtered_aq\
[clc_data_filtered_aq["Count"] >= 2]

# discarded variants 3 (keep discared data)
clc_data_discarded_3 = clc_data_filtered_aq\
[clc_data_filtered_aq["Count"] < 2]

#-------------------------------- 
# Filter data with Frequency >= 5
#--------------------------------
clc_data_filtered_frq = clc_data_filtered_c\
[clc_data_filtered_c["Frequency"] >= 5]

# discarded variants 4 (keep discared data)
clc_data_discarded_4 = clc_data_filtered_c\
[clc_data_filtered_c["Frequency"] < 5]

#--------------------------------------------- 
# Filter data with Forward/reverse balance > 0
#---------------------------------------------
clc_data_filtered_frb = clc_data_filtered_frq\
[clc_data_filtered_frq["Forward/reverse balance"] > 0]

# discarded variants 5 (keep discared data)
clc_data_discarded_5 = clc_data_filtered_frq\
[clc_data_filtered_frq["Forward/reverse balance"] <= 0]

#------------------------------------------------------------ 
# Filter data with Read position test probability >= 0.000001
#------------------------------------------------------------
clc_data_filtered_rptp = clc_data_filtered_frb\
[clc_data_filtered_frb["Read position test probability"] >= 0.000001]

# discarded variants 6 (keep discared data)
clc_data_discarded_6 = clc_data_filtered_frb\
[clc_data_filtered_frb["Read position test probability"] < 0.000001]

#------------------------------------------------------------- 
# Filter data with Read direction test probability >= 0.000001
#-------------------------------------------------------------
clc_data_filtered = clc_data_filtered_rptp\
[clc_data_filtered_rptp["Read direction test probability"] >= 0.000001]

# discarded variants 6 (keep discared data)
clc_data_discarded_7 = clc_data_filtered_rptp\
[clc_data_filtered_rptp["Read direction test probability"] < 0.000001]

#-------------------------
# filter non-synonymous???
#-------------------------

#---------------------------------------
# Generate xlsx file of removed_variants
#---------------------------------------
discarded_clc_data = [clc_data_discarded_1, clc_data_discarded_2, \
                      clc_data_discarded_3, clc_data_discarded_4, \
                      clc_data_discarded_5, clc_data_discarded_6, \
                      clc_data_discarded_7]

# concatenate discarded/removed clc_data variants
removed_variants = pd.concat(discarded_clc_data)

# write file
removed_variants.to_excel(args.removed_variants, \
                            index = False, \
                             engine= None) 

#--------------------------------------------------
# Reindex clc data for further processing necessary
#-------------------------------------------------- 
clc_data_filtered = clc_data_filtered.reset_index()

#------------------------------------------------------------------------
# "Explode" column "coding region change"
# "Explode" necessary since each row contains several transcript:HGVS values
# Result: presentation of one transcript:HGVS per row; a variant shown in a 
# variety of rows derived from the total transcript number 
#------------------------------------------------------------------------           
clc_data_filtered = clc_data_filtered.rename(columns=\
           {"Coding region change": "Coding_region_change"})

clc_data_filtered = (
    clc_data_filtered.assign(Coding_region_change=clc_data_filtered\
                                  ["Coding_region_change"].str.split("; "))
        .explode("Coding_region_change")
        .reset_index(drop=True)
)

#---------------------------------------------------------------------------
# Get NM/XM and HGVSc derived from transcript:HGVS in "Coding_region_change"
# NM/XM and HGVSc are saved in new columns ("NM_v" and "HGVSc)
#---------------------------------------------------------------------------    
for index, NM_tr in enumerate(clc_data_filtered["Coding_region_change"]):
    
    if pd.isna(NM_tr):
        clc_data_filtered.loc[index, "NM_v"] = clc_data_filtered.loc\
                                               [index, "NM_v"]
        clc_data_filtered.loc[index, "HGVSc"] = clc_data_filtered.loc\
                                                [index, "HGVSc"]
        
    else:
        clc_data_filtered.loc[index, "NM_v"] = NM_tr.split(":")[0]
        clc_data_filtered.loc[index, "HGVSc"] = NM_tr.split(":")[1]
        
#------------------------------
# Drop Na values from DataFrame (necessary?)
#------------------------------
clc_data_filtered_dropNa = clc_data_filtered[clc_data_filtered\
                          ["NM_v"].str.contains("nan") == False]

#------------------------------------------
# Load PANCANCER RefSeq transcripts to list (RefSeq check with Natalie und Anna-Lena!!!)
#------------------------------------------
# transcript_list = ".xlsx"
transcript_list = args.transcripts
RefSeq_NM = pd.read_excel(transcript_list)
RefSeq_NM_lst = RefSeq_NM["NM_RefSeq_final"].values.tolist()

#-------------------------------
# Check transcript input for " "
#-------------------------------
for RefSeq_idx in range(len(RefSeq_NM)):
    if ' ' in RefSeq_NM.loc[RefSeq_idx, "NM_RefSeq_final"]:
        raise ValueError('Space in transcript name! Please correct!')
  
#------------------------------------------------------------------------------
# Reset indices; necessary for further processing; since index could be not
# in order due to step "Drop Na values from DataFrame" (see above) (necessary?)
#------------------------------------------------------------------------------
clc_data_filtered_dropNa = clc_data_filtered_dropNa.reset_index()

#----------------------------
# convert to str (necessary?)
#----------------------------
clc_data_filtered_dropNa["NM_v"] = clc_data_filtered_dropNa["NM_v"].astype(str)

#-----------------------------------------------------
# Get index of rows presented in PANCANCER RefSeq list
#-----------------------------------------------------
NM_idx = []
for i in range(len(clc_data_filtered_dropNa["NM_v"])):
    if clc_data_filtered_dropNa["NM_v"][i].split(".")[0] in RefSeq_NM_lst:
        NM_idx.append(i)

#-----------------------------------------
# Filter variants accoriing to NM_idx list
# Store result in new variable
#-----------------------------------------
pre_final_data = clc_data_filtered_dropNa.loc[NM_idx, :]

#------------------------------------------------------------
# Convert NM + version in NM only for subsequent merging step
# Store NM only in column NM_merge
#------------------------------------------------------------
index_variants = pre_final_data.index.tolist()
for i in range(len(pre_final_data["NM_v"])):
    tmp_NM_name = pre_final_data["NM_v"][index_variants[i]].split(".")
    pre_final_data.loc[[index_variants[i]], "NM_merge"] =  tmp_NM_name[0]

#--------------------------------------------
# Generate clinvar Link and add to data 
# Add column "Clinvar_Link" and generate link
#--------------------------------------------
pre_final_data["Clinvar_Link"] = ""
rs = pre_final_data["name dbsnp_v151_ensembl_hg38_no_alt_analysis_set"]
rs_idx = rs.index.tolist()
link_lst = []

#-----------------------
# Get clinvar https list
#-----------------------
pre_final_data["Clinvar_Link"] = fp.clinvar_link_list(link_lst,\
                                 rs_idx, rs)

# Log information
print("Process information:")
print("--> Processing CLC_PanCancer data: successful!")    
    
#-------------------------------------------------
# Process VEP data (obtained via ensembl-vep tool)
#-------------------------------------------------
# vep_file = ".txt"
vep_file = args.vep
VEP_data = pd.read_csv(vep_file, delimiter="\t")  

#-----------------------------------------------------------------                    
# Add column "Chromosome" and "Position"
# Get chromosome and position from column "Location" and NM values
# Result: Ready VEP data for merging with CLC_PAN_data
#-----------------------------------------------------------------     
VEP_data = fp.adjust_chr_pos_NM_VEP(VEP_data)
# ok bis hier 06.02.2024
#------------------------------------
# Merge CLC_PAN_data and VEP data dfs
#------------------------------------

#pre_final_data["Chromosome"] = pre_final_data["Chromosome"] ???
#VEP_data["Chromosome"] = VEP_data["Chromosome"] ???

#----------------------------------------------
# Prepare CLC_PAN_data and VEP data for merging
#----------------------------------------------
VEP_data["Chromosome"] = VEP_data["Chromosome"].astype(str)
pre_final_data["Chromosome"] = pre_final_data["Chromosome"].astype(str)

VEP_data["Position"] = VEP_data["Position"].astype(int)
pre_final_data["Position"] = pre_final_data["Position"].astype(int)
                 
VEP_data = VEP_data.rename\
                        (columns={"Feature": "NM_merge"})
                        
VEP_data["NM_merge"] = VEP_data["NM_merge"].astype(str)
pre_final_data["NM_merge"] = pre_final_data["NM_merge"].astype(str)

print("--> Processing VEP_Ensembl data: successful!")  

#------------------------------------------------------------------------------
# Merge CLC_PAN_data with VEP via columns "Chromosome", "Position", "NM_merge"
# and add aditional columns from vep
#------------------------------------------------------------------------------
merged = pd.merge(pre_final_data, VEP_data[["Chromosome", "Position",\
                 "NM_merge", "HGVSc", "HGVSp", "SYMBOL", "AF", "MAX_AF", \
                 "gnomADe_AF", "gnomADg_AF", "SIFT", "PolyPhen", "CLIN_SIG", \
                 "PUBMED" ]], \
                 on = ["Chromosome", "Position", "NM_merge"], how = "left")
 
#------------------------------
# Get HGVSp nomenclature RefSeq
# Get indices necessary ???
#------------------------------
index_variants_merged = merged.index.tolist()
for i in range(len(merged["HGVSp"])):
    if pd.isna(merged["HGVSp"][index_variants_merged[i]]):
        merged.loc[[index_variants_merged[i]], "HGVS_PROTEIN"] = merged\
            ["HGVSp"][index_variants_merged[i]]
    else:
        tmp_p_name = merged["HGVSp"][index_variants_merged[i]].split(":")
        if len(tmp_p_name) >= 2:
            merged.loc[[index_variants_merged[i]], "HGVS_PROTEIN"] =  tmp_p_name[1]

#-----------------------------------------------------------   
# Merge/join internal variantDB (variantDBi) PAN_CANCER_DATA
# Change to current "Variantenliste" if needed
#------------------------------------------------------------
# lst_var = ".xlsx"
# variantDBi = pd.read_excel(lst_var)
variantDBi = pd.read_excel(args.variant_DBi)

merged = pd.merge(merged,\
                  variantDBi,\
                  left_on = ["name dbsnp_v151_ensembl_hg38_no_alt_analysis_set"],\
                  right_on = ["name dbsnp_v151_ensembl_hg38_no_alt_analysis_set"],\
                  how = "left")

#--------------------------------
# Get comprehensive output format
#--------------------------------
processed_data_final = merged[["Chromosome", "Position", "End Position", \
                       "Reference", "Allele", "Count", "Coverage", \
                        "Frequency", "QUAL", "Forward/reverse balance", \
                        "Average quality","Read position test probability", \
                        "Read direction test probability", "BaseQRankSum",
                        "Homopolymer", "Homopolymer length", \
                        "Count (singleton UMI)", "Count (big UMI)", \
                        "Proportion (singleton UMIs)", "SYMBOL", \
                        "NM_v", "NM_merge", "HGVSc_x", "HGVS_PROTEIN", "Exon Number", \
                        "func dbsnp_v151_ensembl_hg38_no_alt_analysis_set", \
                        "name dbsnp_v151_ensembl_hg38_no_alt_analysis_set", \
                        "CLNSIG clinvar_20220730_hg38_no_alt_analysis_set", \
                        "CLIN_SIG", \
                        "CLNREVSTAT clinvar_20220730_hg38_no_alt_analysis_set", \
                        "AF_EXAC clinvar_20220730_hg38_no_alt_analysis_set", \
                        "AF_TGP clinvar_20220730_hg38_no_alt_analysis_set", \
                        "AF", "MAX_AF", "gnomADe_AF", "gnomADg_AF", "SIFT", \
                        "PolyPhen", "PUBMED", "Wertung", "Clinvar_Link"]]
   
#---------------------------
# Round AF to max 2 decimals 
#---------------------------
processed_data_final.loc[:,"Frequency"]  = processed_data_final\
                                   ["Frequency"].apply(lambda x: round(x,2))

#---------------                                   
# Sort Frequency
#---------------     
final_variants = processed_data_final.sort_values(by=["Frequency"], ascending=False)  

#--------------------
# Identify duplicates
#--------------------  

# Define dict for agg function
agg_dict = {"Chromosome": "first", \
            "Position": "first", \
            "End Position": "first", \
            "Reference": "first", \
            "Allele": "first", \
            "Count": "first", \
            "Coverage": "first", \
            "Frequency": "first", \
            "QUAL": "first", \
            "Forward/reverse balance": "first", \
            "Average quality": "first", \
            "Read position test probability": "first", \
            "Read direction test probability": "first", \
            "BaseQRankSum": "first", \
            "Homopolymer": "first", \
            "Homopolymer length": "first", \
            "Count (singleton UMI)": "first", \
            "Count (big UMI)": "first", \
            "Proportion (singleton UMIs)": "first", \
            "SYMBOL": "first", \
            "NM_v": "first", \
            "NM_merge": "first", \
            "HGVSc_x": "first", \
            "HGVS_PROTEIN": "first", \
            "Exon Number": "first", \
            "func dbsnp_v151_ensembl_hg38_no_alt_analysis_set": "first", \
            "name dbsnp_v151_ensembl_hg38_no_alt_analysis_set": "first", \
            "CLNSIG clinvar_20220730_hg38_no_alt_analysis_set": "first", \
            "CLIN_SIG": "first", \
            "CLNREVSTAT clinvar_20220730_hg38_no_alt_analysis_set": "first", \
            "AF_EXAC clinvar_20220730_hg38_no_alt_analysis_set": "first", \
            "AF_TGP clinvar_20220730_hg38_no_alt_analysis_set": "first", \
            "AF": "first", \
            "MAX_AF": ",".join, \
            "gnomADe_AF": ",".join, \
            "gnomADg_AF": ",".join, \
            "SIFT": "first", \
            "PolyPhen": "first", \
            "PUBMED": "first", \
            "Wertung": "first", \
            "Clinvar_Link": "first"}
    
# list for discarded duplicates dfs
discarded_duplicates = []

# Get string presenting Chromosome_Position_End Postion 
final_variants["check_duplicated"] = ""
col_combine = ["Chromosome", "Position", "End Position"]
final_variants["check_duplicated"] = final_variants[col_combine].\
                                     apply(lambda row: '_'.\
                                     join(row.values.astype(str)), axis=1)
                                                    
# Get duplicate rows
duplicates_rows = final_variants["check_duplicated"][final_variants\
                               ["check_duplicated"].duplicated() == True].tolist()

if duplicates_rows != []:
    
    for i in range(len(duplicates_rows)):
        duplicates_index = final_variants[final_variants\
                 ["check_duplicated"] == duplicates_rows[i]].index.tolist()
            
        if len(duplicates_index) == 2:
        
        # print(duplicates_index)
    
            combine_duplicates = final_variants.loc[duplicates_index, :]
            
            discarded_duplicates.append(combine_duplicates)
        
            tmp_cd = combine_duplicates.groupby\
                             ("check_duplicated", as_index=False).agg(agg_dict)
        
            final_variants.loc[duplicates_index[0], :] = tmp_cd.loc[0, :]
            
            final_variants.drop(index=duplicates_index[1], inplace=True)

    # concatenate discarded/removed duplicated variants
    discarded_duplicates_df = pd.concat(discarded_duplicates)

    # write file
    discarded_duplicates_df.to_excel(args.removed_duplicates, \
                            index = False, \
                             engine= None) 

else:
    
    write_data = ["No duplicates"]
    discarded_duplicates_df = pd.DataFrame(write_data, columns=["No_Duplicates"])
    
    # write file
    discarded_duplicates_df.to_excel(args.removed_duplicates, \
                            index = False, \
                             engine= None) 
            
# Log information
print("--> Combining and formatting CLC and VEP data: successful!")

#----------
# Save file
#----------
final_variants.to_excel(args.outfile, \
                            index = False, \
                             engine= None) 

# Log information
print("--> Writing data to xlsx output files: successful!")
# Getting current date and time for log information
date_time_now = datetime.now()
# dd/mm/YY H:M:S
dt_string = date_time_now.strftime("%d/%m/%Y %H:%M:%S")
print("End:", dt_string)
 
# PLugins
# finale Dokumentation und Validierung





































