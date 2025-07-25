#!/usr/bin/env python

import pandas as pd
import argparse
import json as js
import functions_pan_etl.MVdataset_generator_utils as etl

# using argparse for positinal arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--sample_id", type=str)
parser.add_argument("-x", "--xlsx_path", type=str)
parser.add_argument("-f", "--fastp_json", type=str)
parser.add_argument("-a", "--fq_256_json", type=str)
parser.add_argument("-b", "--fq_bytes_json", type=str)
parser.add_argument("-c", "--bam_json", type=str)
parser.add_argument("-d", "--vcf_256_json", type=str)
parser.add_argument("-e", "--vcf_bytes_json", type=str)
parser.add_argument("-g", "--bed_256_json", type=str)
parser.add_argument("-j", "--bed_bytes_json", type=str)
args = parser.parse_args()
      
# make pandas excel
filepath = args.xlsx_path # xlsx_path !!!
excel_file = pd.ExcelFile(filepath, engine="openpyxl")
    
# get smallVariants (SNVs) as dict records # add coverage, frequency, klinische relevanz,
# Interpretation Bewertung?
if "Auswertung" in excel_file.sheet_names:               
    idx = excel_file.sheet_names.index("Auswertung")
    smallVariatns_records = etl.pan_snvdata_to_dicts\
                            (filepath,excel_file,args.sample_id,idx)
        
elif "auswertung" in excel_file.sheet_names:
    idx = excel_file.sheet_names.index("auswertung")
    smallVariatns_records = etl.pan_snvdata_to_dicts\
                            (filepath,excel_file,args.sample_id,idx)
                                
# get id for check, entity, cellularity, MSI, TMB (more flexible)
if "Final" in excel_file.sheet_names:               
    idx1 = excel_file.sheet_names.index("Final")
    pan_final_page_dict = etl.pan_final_page_to_dict\
                          (filepath,excel_file,idx1)
            
elif "final" in excel_file.sheet_names:
    idx1 = excel_file.sheet_names.index("final")
    pan_final_page_dict = etl.pan_final_page_to_dict\
                          (filepath,excel_file,idx1)
                              
# write to mv_data - use
Oncology_Molecular_Report = { # pancancer data MV + adds
"smallVariatns" : smallVariatns_records,
"copyNumberVariants" : [{
    "identifier" : "not_available",
    "genomicSource" : "not_available",
    "cnvType" : "not_available",
    "gene.display" : "not_available",
    "gene.code" : "not_available",
    "chromosome" : "not_available",
    "endPosition" : "not_available",
    "startPosition" : "not_available",
    "localization" : "not_available"
    }],
"structualVariants" : [{
    "identifier" : "not_available",
    "genomicSource" : "not_available",
    "geneA.display" : "not_available",
    "geneA.code" : "not_available",
    "geneB.diplay" : "not_available",
    "geneB.code" : "not_available",
    "structureType" : "not_available",
    "description" : "not_available"
     }],
"expressionVariants" : [{
    "identifier" : "not_available",
    "gene.display" : "not_available",
    "gene.code" : "not_available",
    "expressionType" : "not_available",
    "reference" : "not_available"
    }],
"complexBiomarkers": [{
    "ploidy" : "not_available",
    "tmb" : pan_final_page_dict["tmb"]["tmb_status"],
    "tmb_mutations_Mb" : pan_final_page_dict["tmb"]["tmb_mutations_Mb"],
    "msi_status" : pan_final_page_dict["msi"]["msi_status"],
    "msiAnzahlMarkerStabil" : pan_final_page_dict["msi"]["msiAnzahlMarkerStabil"],
    "msiAnzahlMarkerInstabil" : pan_final_page_dict["msi"]["msiAnzahlMarkerInstabil"],
    "msiAnzahlMarkerfail" : pan_final_page_dict["msi"]["msiAnzahlMarkerfail"],              
    "hrdHigh" : "not_available",
    "lstHigh" : "not_available",
    "tailHigh" : "not_available"
    }],
"sbsSignatures" : [{
    "sbsVersion" : "not_available",
    "sbsSignatures" : "not_available",
    "sbbsSignaturesPresent" : "not_available"
    }] 
}

# add cellularity
etl.submission_grz["labData"][0]\
                  ["tumorCellCount"][0]\
                  ["count"] = float(pan_final_page_dict["cellularity"]) * 100
                      
# add fastp qc info to subKDK - "percentBasesAboveQualityThreshold"
with open(args.fastp_json, "r") as fastp_json:
    fastp_qc = js.load(fastp_json)
    
reads_q30_rate = round(fastp_qc["summary"]["before_filtering"]["q30_rate"],2) * 100

etl.submission_grz["sequenceData"]\
                  ["percentBasesAboveQualityThreshold"] = reads_q30_rate 
                  
# add bamfile info to subKDK
with open(args.bam_json, "r") as bam_json:
    bam_qc = js.load(bam_json)
   
etl.submission_grz["sequenceData"]\
                  ["meanDepthOfCoverage"] = bam_qc["bam_qc"][0]["mean_cov"]

etl.submission_grz["sequenceData"]\
                  ["minCoverage"] = bam_qc["bam_qc"][0]["min_cov"]

etl.submission_grz["sequenceData"]\
                  ["targetedRegionsAboveMinCoverage"] = bam_qc["bam_qc"][0]\
                                                        ["target_above_mincov"]
        
#add fastq info to submission_grz
with open(args.fq_256_json, "r") as fq_256:
    fq_sha256sum = js.load(fq_256)
    
with open(args.fq_bytes_json, "r") as fq_bytes:
    fq_bytesize = js.load(fq_bytes)
    
# concate/merge json info to dict
fq_sha256sum["fastq_checksums"][0]["R1"]["fileSizeInBytes"] = fq_bytesize\
                                                              ["fastq_bytesizes"]\
                                                              [0]["R1"]\
                                                              ["fileByteSize"]

fq_sha256sum["fastq_checksums"][1]["R2"]["fileSizeInBytes"] = fq_bytesize\
                                                              ["fastq_bytesizes"]\
                                                              [1]["R2"]\
                                                              ["fileByteSize"]

files_dict_cal_concat = {
                        list(fq_sha256sum["fastq_checksums"]\
                        [0].keys())[0]: fq_sha256sum["fastq_checksums"][0]["R1"],
                        list(fq_sha256sum["fastq_checksums"]\
                        [1].keys())[0]: fq_sha256sum["fastq_checksums"][1]["R2"]                       
                     }
                                                        
indices_reads = [0,1]
if len(etl.submission_grz["files"]) == 4:
    # index 0 for read1
    # index 1 for read2
    for index_read in indices_reads:
        print(index_read)
        key, value = list(files_dict_cal_concat.items())[index_read]
          
        etl.submission_grz["files"]\
                          [index_read]["filePath"] = files_dict_cal_concat[key]["file"]
        etl.submission_grz["files"]\
                          [index_read]["fileType"] = "fastq"
        etl.submission_grz["files"]\
                          [index_read]["checksumType"] = "sha256"
        etl.submission_grz["files"]\
                          [index_read]["fileChecksum"] = files_dict_cal_concat[key]["fileChecksum"]
        etl.submission_grz["files"]\
                          [index_read]["fileSizeInBytes"] = files_dict_cal_concat[key]["fileSizeInBytes"]
        etl.submission_grz["files"]\
                          [index_read]["readOrder"] = key
                                 
#add vcf info to submission_grz
with open(args.vcf_256_json, "r") as vcf_256:
    vcf_checksum = js.load(vcf_256)
    
with open(args.vcf_bytes_json, "r") as vcf_bytes:
    vcf_bytesize = js.load(vcf_bytes)
    
if vcf_checksum["vcf_checksum"][0]["file"] == vcf_bytesize["vcf_bytesize"][0]["file"]:
    vcf_name = vcf_checksum["vcf_checksum"][0]["file"]
        
vcf_sha256sum = vcf_checksum["vcf_checksum"][0]["fileChecksum"]

vcf_size = vcf_bytesize["vcf_bytesize"][0]["fileByteSize"]
    
index_vcf = 2
if len(etl.submission_grz["files"]) == 4:
    # index 2 for vcf
    etl.submission_grz["files"]\
                      [index_vcf]["filePath"] = vcf_name
    etl.submission_grz["files"]\
                      [index_vcf]["fileType"] = "vcf"
    etl.submission_grz["files"]\
                      [index_vcf]["checksumType"] = "sha256"
    etl.submission_grz["files"]\
                      [index_vcf]["fileChecksum"] = vcf_sha256sum
    etl.submission_grz["files"]\
                      [index_vcf]["fileSizeInBytes"] =  vcf_size
                          
#add bed info to submission_grz
with open(args.bed_256_json, "r") as bed_256:
    bed_checksum = js.load(bed_256)
    
with open(args.bed_bytes_json, "r") as bed_bytes:
    bed_bytesize = js.load(bed_bytes)
    
if bed_checksum["bedfile_checksum"][0]["file"] == bed_bytesize["bed_bytesize"]\
                                                           [0]["file"]:
    bed_name = bed_checksum["bedfile_checksum"][0]["file"]
        
bed_sha256sum = bed_checksum["bedfile_checksum"][0]["fileChecksum"]

bed_size = bed_bytesize["bed_bytesize"][0]["fileByteSize"]
 
index_bed = 3
if len(etl.submission_grz["files"]) == 4:

    etl.submission_grz["files"]\
                      [index_bed]["filePath"] = bed_name
    etl.submission_grz["files"]\
                      [index_bed]["fileType"] = "bed"
    etl.submission_grz["files"]\
                      [index_bed]["checksumType"] = "sha256"
    etl.submission_grz["files"]\
                      [index_bed]["fileChecksum"] = bed_sha256sum
    etl.submission_grz["files"]\
                      [index_bed]["fileSizeInBytes"] =  bed_size
                                      
#if patient_id == patient_id_check:
pancancer_json = {
       "analysis" : "pancancer",
       "pathoProId_or_nexusId" : "not_yet_available",
       "orbisId_or_orbis_caseId" :  "not_yet_available",
       "network" : "Modellvorhaben_Genomsequezierung",
       "diseaseType" : "oncological",
       "entity" : pan_final_page_dict["entity"],
       "molecularReportod" : Oncology_Molecular_Report,
       "submission_grz" : etl.submission_grz
        }
    
with open(args.sample_id + "_submit.json", "w") as file:
    js.dump(
    pancancer_json, 
    file, 
    indent=4, 
    sort_keys=False
)
        
