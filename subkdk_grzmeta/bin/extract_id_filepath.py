#!/usr/bin/env python

import os
import pandas as pd
import argparse
import re
from pathlib import Path
import glob

# using argparse for positinal arguments
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--target_dir_mvpan", type=str)
parser.add_argument("-a", "--nextseq_data_dir", type=str)
parser.add_argument("-c", "--novaseq_data_dir", type=str)
parser.add_argument("-p", "--pan_ie_dir", type=str)
parser.add_argument("-x", "--id_xlsx_paths", type=str)
parser.add_argument("-f", "--id_fastq_paths", type=str)
parser.add_argument("-b", "--id_bam_path", type=str)
parser.add_argument("-v", "--id_vcf_path", type=str)
args = parser.parse_args()

# ETL path: Pancancer MV
# target dirs Pan
#parent_directory = os.path.abspath('.')
parent_directory = os.path.abspath(args.target_dir_mvpan)
target_dirs = [f for f in os.listdir(parent_directory) if re.match(r'[0-9]', f)]

# dirs done
dirs_done = []
for dirs in target_dirs:
       if str(dirs).endswith("_fertig"):
           dirs_done.append(parent_directory+"/"+dirs)

# get all Pan DONE files from target dir(s)
files_DONE_pan = []
for date_dirs_done in dirs_done:
    print(date_dirs_done)
    for done in os.listdir(date_dirs_done):
        print(done)
        if str(done).endswith("_fertig"):
            check_is_dir = date_dirs_done  + "/" +done
            if not os.path.isdir(check_is_dir):
                raise ValueError(f"'{check_is_dir}' is not a directory.")

            paths = glob.glob(date_dirs_done  + "/" + done +  "/*/FINAL_OUTPUT/*_final_processed_P.xlsx")
            for path in paths:
                files_DONE_pan.append(str(path))

# make id_path_df for json generator
size_of_df = len(files_DONE_pan)
id_path_df = pd.DataFrame("", index=list(range(0,size_of_df)),
                          columns = ["patient_id", "xlsx_path", "fq_R1_path",
                                     "fq_R2_path", "bam_path", #"bai_path",
                                     "vcf_path"])

for i, pancancer_file in enumerate(files_DONE_pan):
    print(i, pancancer_file)

    # get case id and set xlsx_path
    basename = os.path.basename(pancancer_file)
    patient_id = basename.split("_")[0] # get check !!!
    id_path_df.loc[i,"patient_id"] = str(patient_id)
    id_path_df.loc[i,"xlsx_path"] = Path(pancancer_file)

    # find fastq.gz files of patient - function
    raw_data_dir_NextSeq = args.nextseq_data_dir
    raw_data_dir_NovaSeq = args.novaseq_data_dir

    fastq_read1_inNextSeq = glob.glob(raw_data_dir_NextSeq +"/*/"+patient_id+"*_1.fq.gz")
    fastq_read2_inNextSeq = glob.glob(raw_data_dir_NextSeq +"/*/"+patient_id+"*_2.fq.gz")
    fastq_read1_inNovaSeq = glob.glob(raw_data_dir_NovaSeq+"/*/"+patient_id+"*_1.fq.gz")
    fastq_read2_inNovaSeq= glob.glob(raw_data_dir_NovaSeq+"/*/"+patient_id+"*_2.fq.gz")

    if fastq_read1_inNextSeq != [] and fastq_read1_inNextSeq != []:
        path_fastq_read1 = fastq_read1_inNextSeq
        path_fastq_read2 = fastq_read2_inNextSeq

    elif fastq_read1_inNovaSeq != [] and fastq_read1_inNovaSeq != []:
        path_fastq_read1 = fastq_read1_inNovaSeq
        path_fastq_read2 = fastq_read2_inNovaSeq

    id_path_df.loc[i,"fq_R1_path"] = Path(path_fastq_read1[0])
    id_path_df.loc[i,"fq_R2_path"] = Path(path_fastq_read2[0])

    # find bamfiles
    bam_file_dir = args.pan_ie_dir
    if glob.glob(bam_file_dir +"/**/"+patient_id+".bam") != []:
        bam_file_path = glob.glob(bam_file_dir +"/**/"+patient_id+".bam")
    elif glob.glob(bam_file_dir +"/**/**/"+patient_id+".bam") != []:
        bam_file_path = glob.glob(bam_file_dir +"/**/**/"+patient_id+".bam")

    id_path_df.loc[i,"bam_path"] = Path(bam_file_path[0])
    #id_path_df.loc[i,"bai_path"] = Path(bam_file_path[0]+".bai")

    # find vcf file PanCancer
    # check that only one file exits! also for fastq / all data
    vcf_files_dir = args.pan_ie_dir
    if glob.glob( vcf_files_dir  +"/**/"+patient_id+".vcf") != []:
        vcf_file_path = glob.glob( vcf_files_dir  +"/**/"+patient_id+".vcf")
    elif glob.glob( vcf_files_dir  +"/**/**/"+patient_id+".vcf") != []:
        vcf_file_path = glob.glob( vcf_files_dir  +"/**/**/"+patient_id+".vcf")

    id_path_df.loc[i,"vcf_path"] = Path(vcf_file_path[0])

# make samplesheets
# for qc_fastp and qc_samtools_depth
#fastp
qc_fastp_df = id_path_df[["patient_id", "fq_R1_path","fq_R2_path"]]
qc_fastp_df.to_csv(args.id_fastq_paths, sep=",",index=False)
#samtools depth
#qc_samtools_depth_df = id_path_df[["patient_id", "bam_path", "bai_path"]]
qc_samtools_depth_df = id_path_df[["patient_id", "bam_path"]]
qc_samtools_depth_df.to_csv(args.id_bam_path, sep=",",index=False)
# for vcf calculations
vcf_df = id_path_df[["patient_id", "vcf_path"]]
vcf_df.to_csv(args.id_vcf_path, sep=",",index=False)
# for json generator
id_xlsx_df = id_path_df[["patient_id", "xlsx_path"]]
id_xlsx_df.to_csv(args.id_xlsx_paths, sep=",",index=False)
