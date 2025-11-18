#!/usr/bin/env python

import os
import pandas as pd
import argparse
import re
from pathlib import Path
import glob

# using argparse for positinal arguments
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--target_dir_mvwgs", type=str)
#parser.add_argument("-a", "--nextseq_data_dir", type=str)
parser.add_argument("-c", "--novaseq_data_dir", type=str)
parser.add_argument("-p", "--nxf_outputdir", type=str)
parser.add_argument("-x", "--id_xlsx_paths", type=str)
parser.add_argument("-f", "--id_fastq_paths", type=str)
parser.add_argument("-b", "--id_bam_path", type=str)
parser.add_argument("-v", "--id_vcf_path", type=str)
parser.add_argument("-i", "--id_patient", type=str)
args = parser.parse_args()

# ETL path: WGS MV
#parent_directory = os.path.abspath('.')
parent_directory = os.path.abspath(args.target_dir_mvwgs)
#target_dirs = [f for f in os.listdir(parent_directory) if re.match(r'[0-9]', f)]
files_DONE_wgs = glob.glob(parent_directory + "/" + "*_fertig.xlsx", recursive=True)

# make id_path_df for json generator
size_of_df = len(files_DONE_wgs)
id_path_df = pd.DataFrame("", index=list(range(0,size_of_df)),
                          columns = ["patient_id", 
                                     "xlsx_path", 
                                     "fq_R1_path_n",
                                     "fq_R2_path_n", 
                                     "fq_R1_path_t", 
                                     "fq_R2_path_t",
                                     "bam_path_n", 
                                     "bai_path_n",
                                     "bam_path_t", 
                                     "bai_path_t",
                                     "vcf_path"])


for i, wgs_file in enumerate(files_DONE_wgs):
    print(i, wgs_file)
    
    # get case id and set xlsx_path
    basename = os.path.basename(wgs_file)
    patient_id = basename.split("_")[0] # get check !!!
    # check
    if basename.count(patient_id) == 2:
        id_path_df.loc[i,"patient_id"] = str(patient_id)
    else:
        raise ValueError(f"Please check the following patient '{patient_id}'")
        
    id_path_df.loc[i,"xlsx_path"] = Path(wgs_file)
        
    # find fastq.gz files of patient with correct tissue
    #raw_data_dir_NextSeq = args.nextseq_data_dir
    raw_data_dir_NovaSeq = args.novaseq_data_dir
    #raw_data_dir_NextSeq = nextseq_data_dir
    #raw_data_dir_NovaSeq = novaseq_data_dir

    # own name schema - normal
    #fastq_read1_inNextSeq_n = glob.glob(raw_data_dir_NextSeq +"/**/"+patient_id+"*N*_1.fq.gz", recursive=True)
    #fastq_read2_inNextSeq_n = glob.glob(raw_data_dir_NextSeq +"/**/"+patient_id+"*N*_2.fq.gz", recursive=True)
    fastq_read1_inNovaSeq_n = glob.glob(raw_data_dir_NovaSeq+"/"+patient_id+"*N*_1.fq.gz", recursive=True)
    fastq_read2_inNovaSeq_n = glob.glob(raw_data_dir_NovaSeq+"/"+patient_id+"*N*_2.fq.gz", recursive=True)

    #if fastq_read1_inNextSeq_n == [] and fastq_read2_inNextSeq_n == [] and\
    if fastq_read1_inNovaSeq_n == [] and fastq_read2_inNovaSeq_n == []:
        # CF name schema - normal
        #fastq_read1_inNextSeq_n = glob.glob(raw_data_dir_NextSeq +"/**/"+patient_id+"*N*_R1_*.fastq.gz", recursive=True)
        #fastq_read2_inNextSeq_n = glob.glob(raw_data_dir_NextSeq +"/**/"+patient_id+"*N*_R2_*.fastq.gz", recursive=True)
        fastq_read1_inNovaSeq_n = glob.glob(raw_data_dir_NovaSeq+"/"+patient_id+"*N*_R1_*.fastq.gz", recursive=True)
        fastq_read2_inNovaSeq_n = glob.glob(raw_data_dir_NovaSeq+"/"+patient_id+"*N*_R2_*.fastq.gz", recursive=True)


    if fastq_read1_inNovaSeq_n == [] and fastq_read2_inNovaSeq_n == []:
        fastq_read1_inNovaSeq_n = glob.glob(raw_data_dir_NovaSeq+"/"+patient_id+"*B*_1.fq.gz", recursive=True)
        fastq_read2_inNovaSeq_n = glob.glob(raw_data_dir_NovaSeq+"/"+patient_id+"*B*_2.fq.gz", recursive=True)

    if fastq_read1_inNovaSeq_n == [] and fastq_read2_inNovaSeq_n == []:
        fastq_read1_inNovaSeq_n = glob.glob(raw_data_dir_NovaSeq+"/"+patient_id+"*B*_R1_*.fastq.gz", recursive=True)
        fastq_read2_inNovaSeq_n = glob.glob(raw_data_dir_NovaSeq+"/"+patient_id+"*B*_R2_*.fastq.gz", recursive=True)

    # own name schema - tumor
    #fastq_read1_inNextSeq_t = glob.glob(raw_data_dir_NextSeq +"/**/"+patient_id+"*T*_1.fq.gz", recursive=True)
    #fastq_read2_inNextSeq_t = glob.glob(raw_data_dir_NextSeq +"/**/"+patient_id+"*T*_2.fq.gz", recursive=True)
    fastq_read1_inNovaSeq_t = glob.glob(raw_data_dir_NovaSeq+"/"+patient_id+"*T*_1.fq.gz", recursive=True)
    fastq_read2_inNovaSeq_t = glob.glob(raw_data_dir_NovaSeq+"/"+patient_id+"*T*_2.fq.gz", recursive=True)
    
    #if fastq_read1_inNextSeq_t == [] and fastq_read2_inNextSeq_t == [] and\
    if fastq_read1_inNovaSeq_t == [] and fastq_read2_inNovaSeq_t == []:
        # CF name schema - tumor
        #fastq_read1_inNextSeq_t = glob.glob(raw_data_dir_NextSeq +"/**/"+patient_id+"*T*_R1_*.fastq.gz", recursive=True)
        #fastq_read2_inNextSeq_t = glob.glob(raw_data_dir_NextSeq +"/**/"+patient_id+"*T*_R2_*.fastq.gz", recursive=True)
        fastq_read1_inNovaSeq_t = glob.glob(raw_data_dir_NovaSeq+"/"+patient_id+"*T*_R1_*.fastq.gz", recursive=True)
        fastq_read2_inNovaSeq_t = glob.glob(raw_data_dir_NovaSeq+"/"+patient_id+"*T*_R2_*.fastq.gz", recursive=True)
    
    if fastq_read1_inNovaSeq_t == [] and fastq_read2_inNovaSeq_t == []:
        fastq_read1_inNovaSeq_t = glob.glob(raw_data_dir_NovaSeq+"/"+patient_id+"*FF*_1.fq.gz", recursive=True)
        fastq_read2_inNovaSeq_t = glob.glob(raw_data_dir_NovaSeq+"/"+patient_id+"*FF*_2.fq.gz", recursive=True)
    if fastq_read1_inNovaSeq_t == [] and fastq_read2_inNovaSeq_t == []:
        fastq_read1_inNovaSeq_t = glob.glob(raw_data_dir_NovaSeq+"/"+patient_id+"*FF*_R1_*.fastq.gz", recursive=True)
        fastq_read2_inNovaSeq_t = glob.glob(raw_data_dir_NovaSeq+"/"+patient_id+"*FF*_R2_*.fastq.gz", recursive=True)
    #if fastq_read1_inNextSeq_n != [] and fastq_read2_inNextSeq_n != [] and\
    #    fastq_read1_inNextSeq_t != [] and fastq_read2_inNextSeq_t != []:
    #    path_fastq_read1_n = fastq_read1_inNextSeq_n
    #    path_fastq_read2_n = fastq_read2_inNextSeq_n
    #    path_fastq_read1_t = fastq_read1_inNextSeq_t
    #    path_fastq_read2_t = fastq_read2_inNextSeq_t
    # elif
    if fastq_read1_inNovaSeq_n != [] and fastq_read2_inNovaSeq_n != [] and\
        fastq_read1_inNovaSeq_t != [] and fastq_read2_inNovaSeq_t != []:
        path_fastq_read1_n = fastq_read1_inNovaSeq_n
        path_fastq_read2_n = fastq_read2_inNovaSeq_n
        path_fastq_read1_t = fastq_read1_inNovaSeq_t
        path_fastq_read2_t = fastq_read2_inNovaSeq_t
        
    id_path_df.loc[i,"fq_R1_path_n"] = Path(path_fastq_read1_n[0])
    id_path_df.loc[i,"fq_R2_path_n"] = Path(path_fastq_read2_n[0])
    id_path_df.loc[i,"fq_R1_path_t"] = Path(path_fastq_read1_t[0])
    id_path_df.loc[i,"fq_R2_path_t"] = Path(path_fastq_read2_t[0])
    
    # find bamfiles
    # get tissue conservation
    pattern = ["_N_BLOOD", "_N_FFPE", "_N_FF", "_T_FFPE", "_T_FF"]
    # create a regex pattern that matches any of the string
    tissue_conservation_pattern = "|".join(map(re.escape,pattern))
    matches_in_filename = re.findall(tissue_conservation_pattern,basename)
    #nxf_dir = r"W:\ngs_pipeline_data\pipeline_user\nextflow_outputdir"
    nxf_dir = args.nxf_outputdir
    
    # normal
    for tc_pattern in matches_in_filename:
        #print(tc_pattern)
        if "N" in tc_pattern:
            if glob.glob(nxf_dir +"/"+patient_id + tc_pattern + "*.bam", 
               recursive=True) != []:
                bam_file_path_n = glob.glob\
                (nxf_dir +"/"+patient_id + tc_pattern + "*.bam", 
                recursive=True)
            #elif glob.glob(nxf_dir +"/**/**/"+patient_id + tc_pattern + "*.bam", 
            #     recursive=True) != []:
            #    bam_file_path_n = glob.\
            #    glob(nxf_dir +"/**/**/"+patient_id + tc_pattern + "*.bam", 
            #    recursive=True)
    # tumor
        elif "T" in tc_pattern:
            if glob.glob(nxf_dir +"/"+patient_id + tc_pattern + "*.bam", 
               recursive=True) != []:
                bam_file_path_t = glob.glob\
                (nxf_dir +"/"+patient_id + tc_pattern + "*.bam", 
                recursive=True)
            #elif glob.glob(nxf_dir +"/**/**/"+patient_id + tc_pattern + "*.bam", 
            #     recursive=True) != []:
            #    bam_file_path_t = glob.\
            #    glob(nxf_dir +"/**/**/"+patient_id + tc_pattern + "*.bam", 
            #    recursive=True)
                
    id_path_df.loc[i,"bam_path_n"] = Path(bam_file_path_n[0])
    id_path_df.loc[i,"bai_path_n"] = Path(bam_file_path_n[0]+".bai")
    id_path_df.loc[i,"bam_path_t"] = Path(bam_file_path_t[0])
    id_path_df.loc[i,"bai_path_t"] = Path(bam_file_path_t[0]+".bai")
    
    # find vcf file wgs
    # check that only one file exits! also for fastq / all data
    vcf_files_dir = nxf_dir
    if glob.glob( vcf_files_dir  +"/"+patient_id+"*.vcf", recursive=True) != []:
        vcf_file_path = glob.glob( vcf_files_dir  +"/"+patient_id+"*.vcf", recursive=True)
    #elif glob.glob( vcf_files_dir  +"/**/**/"+patient_id+"*.vcf", recursive=True) != []:
    #    vcf_file_path = glob.glob( vcf_files_dir  +"/**/**/"+patient_id+"*.vcf", recursive=True)

    # due to tissue conservation more vcf files can exist - find the correct one
    # using matches_in_filename 
    vcf_file_path_final = []
    if len(matches_in_filename) == 2:
        for vcf_files in vcf_file_path:
            if matches_in_filename[0] in vcf_files and \
               matches_in_filename[1] in vcf_files:
               vcf_file_path_final.append(vcf_files)
   
    id_path_df.loc[i,"vcf_path"] = Path(vcf_file_path_final[0])
    
# make samplesheets
# for qc_fastp and qc_samtools_depth
#fastp
qc_fastp_df = id_path_df[["patient_id", "fq_R1_path_n","fq_R2_path_n",
                              "fq_R1_path_t","fq_R2_path_t"]]
qc_fastp_df.to_csv(args.id_fastq_paths, sep=",",index=False)
#samtools depth
#qc_samtools_depth_df = id_path_df[["patient_id", "bam_path", "bai_path"]]
qc_samtools_depth_df = id_path_df[["patient_id", "bam_path_n", "bai_path_n",
                                       "bam_path_t", "bai_path_t"]]
qc_samtools_depth_df.to_csv(args.id_bam_path, sep=",",index=False)
# for vcf calculations
vcf_df = id_path_df[["patient_id", "vcf_path"]]
vcf_df.to_csv(args.id_vcf_path, sep=",",index=False)
# for json generator
id_xlsx_df = id_path_df[["patient_id", "xlsx_path"]]
id_xlsx_df.to_csv(args.id_xlsx_paths, sep=",",index=False)
# for meta data fmrest
id_patient_df = id_path_df["patient_id"]
id_patient_df.to_csv(args.id_patient, sep=",",index=False)
