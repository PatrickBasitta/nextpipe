#!/usr/bin/env python

import pandas as pd
import os.path
import argparse
import json as js
import functions_etl.MVdataset_generator_utils as etl
import functions_etl.global_variables as gv
import requests
#import datetime

# using argparse for positinal arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--sample_id", type=str)
parser.add_argument("-x", "--xlsx_path", type=str)
parser.add_argument("-f", "--fastp_json_normal", type=str)
parser.add_argument("-a", "--fq_sha256_json_normal", type=str)
parser.add_argument("-b", "--fq_bytes_json_normal", type=str)
parser.add_argument("-r", "--fastp_json_tumor", type=str)
parser.add_argument("-s", "--fq_sha256_json_tumor", type=str)
parser.add_argument("-t", "--fq_bytes_json_tumor", type=str)
parser.add_argument("-c", "--bam_json_normal", type=str)
parser.add_argument("-u", "--bam_json_tumor", type=str)
parser.add_argument("-d", "--vcf_sha256_json", type=str)
parser.add_argument("-e", "--vcf_bytes_json", type=str)
#parser.add_argument("-g", "--bed_sha256_json", type=str)
#parser.add_argument("-j", "--bed_bytes_json", type=str)
parser.add_argument("-p", "--patient_data_json", type=str)
parser.add_argument("-q", "--hgnc", type=str)
parser.add_argument("-w", "--targeted_bam_json_normal", type=str)
parser.add_argument("-z", "--targeted_bam_json_tumor", type=str)
args = parser.parse_args()
      
# make pandas excel
filepath = args.xlsx_path # xlsx_path !!!
excel_file = pd.ExcelFile(filepath, engine="openpyxl")
hgnc_file = args.hgnc
hgnc = pd.read_csv(hgnc_file, sep="\t")    
# get smallVariants (SNVs) as dict records # add coverage, frequency, klinische relevanz,
# Interpretation Bewertung?
if "Auswertung" in excel_file.sheet_names:               
    idx = excel_file.sheet_names.index("Auswertung")
    smallVariatns_records = etl.wxs_snvdata_to_dicts\
                            (filepath,excel_file,args.sample_id,idx,hgnc) # sampleid
        
elif "auswertung" in excel_file.sheet_names:
    idx = excel_file.sheet_names.index("auswertung")
    smallVariatns_records = etl.wxs_snvdata_to_dicts\
                            (filepath,excel_file,args.sample_id,idx.hgnc)
                                
# get id for check, entity, cellularity, MSI, TMB (more flexible)
if "Final" in excel_file.sheet_names:               
    idx1 = excel_file.sheet_names.index("Final")
    wgs_final_page_dict = etl.wxs_final_page_to_dict\
                          (filepath,excel_file,idx1)
            
elif "final" in excel_file.sheet_names:
    idx1 = excel_file.sheet_names.index("final")
    wgs_final_page_dict = etl.wxs_final_page_to_dict\
                          (filepath,excel_file,idx1)
                              
# write to mv_data - use
Oncology_Molecular_Report = { # wgs data MV + adds
"smallVariatns" : smallVariatns_records,
"copyNumberVariants" : [{
    "identifier" : "na",
    "genomicSource" : "na",
    "cnvType" : "na",
    "gene.display" : "na",
    "gene.code" : "na",
    "chromosome" : "na",
    "endPosition" : "na",
    "startPosition" : "na",
    "localization" : "na"
    }],
"structualVariants" : [{
    "identifier" : "na",
    "genomicSource" : "na",
    "geneA.display" : "na",
    "geneA.code" : "na",
    "geneB.diplay" : "na",
    "geneB.code" : "na",
    "structureType" : "na",
    "description" : "na"
     }],
"expressionVariants" : [{
    "identifier" : "na",
    "gene.display" : "na",
    "gene.code" : "na",
    "expressionType" : "na",
    "reference" : "na"
    }],
"complexBiomarkers": [{
    "ploidy" : "na",
    "tmb" : wgs_final_page_dict["tmb"]["tmb_status"],
    "tmb_mutations_Mb" : float(wgs_final_page_dict["tmb"]["mutations_Mb"]),
    "tmb_Anzahl_Mutationen_missense" : int(wgs_final_page_dict["tmb"]["Anzahl_Mutationen_missense"]),
    "msi_status" : wgs_final_page_dict["msi"]["msi_status"],
    "msiErgebnis_MSIsensor_pro" : float(wgs_final_page_dict["msi"]["Ergebnis_MSIsensor_pro"]),
    #"HR_deficiency_score_OA" : wgs_final_page_dict["HR_deficiency_score_OA"],
    "hrdHigh" : "na",
    "lstHigh" : "na",
    "tailHigh" : "na"
    }],
"sbsSignatures" : [{
    "sbsVersion" : "na",
    "sbsSignatures" : "na",
    "sbbsSignaturesPresent" : "na"
    }]
}

# add patient data
with open(args.patient_data_json, "r") as patient_data:
    p_data = js.load(patient_data)

# submission_grz
etl.wgs_submission_grz["submission"]["submissionDate"] = "generated upon upload" #datetime.datetime.now().replace(hour=0, minute=0, second=0, microsecond=0).isoformat().split("T")[0]
etl.wgs_submission_grz["submission"]["submissionType"] = gv.submissionType # test modus
etl.wgs_submission_grz["submission"]["tanG"] = ""
etl.wgs_submission_grz["submission"]["localCaseId"] = ""
etl.wgs_submission_grz["submission"]["coverageType"] = gv.coverageType
etl.wgs_submission_grz["submission"]["submitterId"] = gv.submitterId
etl.wgs_submission_grz["submission"]["genomicDataCenterId"] = gv.genomicDataCenterId
etl.wgs_submission_grz["submission"]["clinicalDataNodeId"] = gv.clinicalDataNodeId
etl.wgs_submission_grz["submission"]["diseaseType"] = gv.diseaseType # assumed
etl.wgs_submission_grz["submission"]["genomicStudyType"] = gv.genomicStudyType # assumed
etl.wgs_submission_grz["submission"]["genomicStudySubtype"] = gv.wxs_genomicStudySubtype
etl.wgs_submission_grz["submission"]["labName"] = gv.labName

# add broad consent information
etl.wgs_submission_grz["donors"][0]["donorPseudonym"] = ""
etl.wgs_submission_grz["donors"][0]["gender"] = p_data["gender"]
etl.wgs_submission_grz["donors"][0]["relation"] = gv.wgs_relation
#etl.wgs_submission_grz["donors"][0]["mvConsent"]["presentationDate"] = ""
#etl.wgs_submission_grz["donors"][0]["mvConsent"]["version"] = ""
#etl.wgs_submission_grz["donors"][0]["mvConsent"]["scope"][0]["type"] = ""
#etl.wgs_submission_grz["donors"][0]["mvConsent"]["scope"][0]["date"] = ""
#etl.wgs_submission_grz["donors"][0]["mvConsent"]["scope"][0]["domain"] = ""
#etl.wgs_submission_grz["donors"][0]["researchConsents"][0]["schemaVersion"] = ""
#etl.wgs_submission_grz["donors"][0]["researchConsents"][0]["presentationDate"] = ""
#etl.wgs_submission_grz["donors"][0]["researchConsents"][0]["scope"] = ""

# fastp info
with open(args.fastp_json_normal, "r") as fastp_json_normal:
    fastp_qc_normal = js.load(fastp_json_normal)

with open(args.fastp_json_tumor, "r") as fastp_json_tumor:
    fastp_qc_tumor = js.load(fastp_json_tumor)

fastp_qc_normal_tumor = [fastp_qc_normal, fastp_qc_tumor]

# bam info
with open(args.bam_json_normal, "r") as bam_json_normal:
    bam_qc_normal = js.load(bam_json_normal)
    
with open(args.bam_json_tumor, "r") as bam_json_tumor:
    bam_qc_tumor = js.load(bam_json_tumor)
    
bam_qc_normal_tumor = [bam_qc_normal, bam_qc_tumor]

# bam info targeted
with open(args.targeted_bam_json_normal, "r") as targeted_bam_json_normal:
    targeted_bam_qc_normal = js.load(targeted_bam_json_normal)

with open(args.targeted_bam_json_tumor, "r") as targeted_bam_json_tumor:
    targeted_bam_qc_tumor = js.load(targeted_bam_json_tumor)

targeted_bam_qc_normal_tumor = [targeted_bam_qc_normal, targeted_bam_qc_tumor]

#add fastq info to submission_grz
with open(args.fq_sha256_json_normal, "r") as fq_sha256_json_normal:
    fq_sha256sum_normal = js.load(fq_sha256_json_normal)
    
with open(args.fq_sha256_json_tumor, "r") as fq_sha256_json_tumor:
    fq_sha256sum_tumor = js.load(fq_sha256_json_tumor)
    
fq_sha256sum_normal_tumor = [fq_sha256sum_normal, fq_sha256sum_tumor ]

with open(args.fq_bytes_json_normal, "r") as fq_byte_json_normal:
    fq_bytesize_normal = js.load(fq_byte_json_normal)
    
with open(args.fq_bytes_json_tumor, "r") as fq_byte_json_tumor:
    fq_bytesize_tumor = js.load(fq_byte_json_tumor)
    
fq_bytesize_normal_tumor = [fq_bytesize_normal, fq_bytesize_tumor ]

# index_number_files: len(2 files in normal, 3 files in tumor)
index_num_files = [2,3]
sampleconservation = [wgs_final_page_dict["sampleConservation_N"],wgs_final_page_dict["sampleConservation_T"]]
# index 0 = normal, index 1 = tumor
index_normal_tumor = [0,1]
for i_nt in index_normal_tumor:
  
    # add info to sumbmission_grz - labData
    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                      ["labDataName"] = gv.wgs_labDataName[i_nt] #p_data["pathoProId"]

    # add info tissue ontology
    if sampleconservation[i_nt] == "blood":
        
        etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                          ["tissueOntology"]["name"] = "UBERON"

        etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                          ["tissueOntology"]["version"] = "1.5"

        etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                          ["tissueTypeId"] = "UBERON:0000178"

        etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                          ["tissueTypeName"] = "blood"

        etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                          ["sampleConservation"] = "fresh-tissue"

    else:
        etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                          ["tissueOntology"]["name"] = gv.wes_tissueOntology_name[i_nt]

        etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                          ["tissueOntology"]["version"] = gv.wes_tissueOntology_version[i_nt]

        etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                          ["tissueTypeId"] = ""

        etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                          ["tissueTypeName"] = ""

        etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                          ["sampleConservation"] = sampleconservation[i_nt]

    # add further labData information
    sampledate = [wgs_final_page_dict["sampledate_N"], wgs_final_page_dict["sampledate_T"]]

    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                      ["sampleDate"] = sampledate[i_nt]

    #etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
    #                  ["sampleConservation"] = sampleconservation[i_nt]

    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                      ["sequenceType"] = gv.wgs_sequenceType

    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                      ["sequenceSubtype"] = gv.wgs_sequenceSubtype[i_nt]

    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                      ["fragmentationMethod"] = gv.wgs_fragmentationMethod #known

    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                      ["libraryType"] = gv.wgs_libraryType # get from excel as well

    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                      ["libraryPrepKit"] = gv.wgs_libraryPrepKit

    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                      ["libraryPrepKitManufacturer"] = gv.wgs_libraryPrepKitManufacturer

    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                      ["sequencerModel"] = wgs_final_page_dict["sequencer"] # or NextSeq550?

    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                      ["sequencerManufacturer"] = gv.wgs_sequencerManufacturer # known

    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                      ["kitName"] = wgs_final_page_dict["kit_name"]

    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                      ["kitManufacturer"] = gv.wgs_kitManufacturer

    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                      ["enrichmentKitManufacturer"] = gv.wgs_enrichmentKitManufacturer

    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                      ["enrichmentKitDescription"] = gv.wgs_enrichmentKitDescription

    barcode = [wgs_final_page_dict["barcode_N"], wgs_final_page_dict["barcode_T"]]

    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                      ["barcode"] = barcode[i_nt]

    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]\
                      ["sequencingLayout"] = gv.wgs_sequencingLayout # known

    # add sequenceData info to grz_submission
    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]["sequenceData"]\
                      ["bioinformaticsPipelineName"] = gv.wgs_bioinformaticsPipelineName

    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]["sequenceData"]\
                      ["bioinformaticsPipelineVersion"] = gv.wgs_bioinformaticsPipelineVersion

    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]["sequenceData"]\
                      ["referenceGenome"] = gv.wgs_referenceGenome

    # add fastp qc info to submission_grz - "percentBasesAboveQualityThreshold"
    #with open(args.fastp_json, "r") as fastp_json:
    #    fastp_qc = js.load(fastp_json)
    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]["sequenceData"]\
                      ["percentBasesAboveQualityThreshold"]["minimumQuality"] = gv.wgs_minimumQuality[i_nt]

    reads_q30_rate = round(fastp_qc_normal_tumor[i_nt]["summary"]["before_filtering"]["q30_rate"],2) * 100

    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]["sequenceData"]\
                      ["percentBasesAboveQualityThreshold"]["percent"] = round(float(reads_q30_rate),0)

    # add bamfile info to submission_grz
    #with open(args.bam_json, "r") as bam_json:
    #    bam_qc = js.load(bam_json)

    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]["sequenceData"]\
                      ["meanDepthOfCoverage"] = round(float(bam_qc_normal_tumor[i_nt]["bam_qc"][0]["mean_cov"]),2)

    #etl.wgs_submission_grz["donors"][0]["labData"][i_nt]["sequenceData"]\
    #                  ["minCoverage"] = round(float(bam_qc_normal_tumor[i_nt]["bam_qc"][0]["min_cov"]),0) # minCov of the BfarmQs
    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]["sequenceData"]\
                       ["minCoverage"] = gv.wgs_minCoverage[i_nt]
    
    # must be targeted
    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]["sequenceData"]\
                      ["targetedRegionsAboveMinCoverage"] = float(targeted_bam_qc_normal_tumor[i_nt]["targeted_bam_qc"][0]\
                                                        ["targets_above_mincov"])

    # add info of nonCodingVariants and caller information
    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]["sequenceData"]\
                      ["nonCodingVariants"] = gv.wgs_nonCodingVariants

    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]["sequenceData"]\
                      ["callerUsed"][0]["name"] = gv.wgs_callerUsed_name

    etl.wgs_submission_grz["donors"][0]["labData"][i_nt]["sequenceData"]\
                      ["callerUsed"][0]["version"] = gv.wgs_callerUsed_version

    # the relative path of the grz_cli tool is always here sample_id/files/file
    #rel_path_grz_cli = args.sample_id+"/files/"
    #rel_path_grz_cli = args.sample_id+"/files/"

    #add fastq info to submission_grz
    #with open(args.fq_256_json, "r") as fq_256:
    #    fq_sha256sum = js.load(fq_256)

    #with open(args.fq_bytes_json, "r") as fq_bytes:
    #    fq_bytesize = js.load(fq_bytes)

    # concate/merge json info to dict
    fq_sha256sum_normal_tumor[i_nt]["fastq_checksums"][0]["R1"]["fileSizeInBytes"] = fq_bytesize_normal_tumor[i_nt]\
                                                              ["fastq_bytesizes"]\
                                                              [0]["R1"]\
                                                              ["fileByteSize"]

    fq_sha256sum_normal_tumor[i_nt]["fastq_checksums"][1]["R2"]["fileSizeInBytes"] = fq_bytesize_normal_tumor[i_nt]\
                                                              ["fastq_bytesizes"]\
                                                              [1]["R2"]\
                                                              ["fileByteSize"]

    files_dict_cal_concat = {
                        list(fq_sha256sum_normal_tumor[i_nt]["fastq_checksums"]\
                        [0].keys())[0]: fq_sha256sum_normal_tumor[i_nt]["fastq_checksums"][0]["R1"],
                        list(fq_sha256sum_normal_tumor[i_nt]["fastq_checksums"]\
                        [1].keys())[0]: fq_sha256sum_normal_tumor[i_nt]["fastq_checksums"][1]["R2"]
                        }

    indices_reads = [0,1]
    if len(etl.wgs_submission_grz["donors"][0]["labData"][i_nt]["sequenceData"]["files"]) == index_num_files[i_nt]:
        # index 0 for read1
        # index 1 for read2
        for index_read in indices_reads:
            print(index_read)
            key, value = list(files_dict_cal_concat.items())[index_read]

            etl.wgs_submission_grz["donors"][0]["labData"][i_nt]["sequenceData"]["files"]\
                              [index_read]["filePath"] = files_dict_cal_concat[key]["file"]
            etl.wgs_submission_grz["donors"][0]["labData"][i_nt]["sequenceData"]["files"]\
                              [index_read]["fileType"] = gv.wgs_files_fastq
            etl.wgs_submission_grz["donors"][0]["labData"][i_nt]["sequenceData"]["files"]\
                              [index_read]["checksumType"] = gv.wgs_checksumType
            etl.wgs_submission_grz["donors"][0]["labData"][i_nt]["sequenceData"]["files"]\
                              [index_read]["fileChecksum"] = files_dict_cal_concat[key]["fileChecksum"]
            etl.wgs_submission_grz["donors"][0]["labData"][i_nt]["sequenceData"]["files"]\
                              [index_read]["fileSizeInBytes"] = int(files_dict_cal_concat[key]["fileSizeInBytes"])
            etl.wgs_submission_grz["donors"][0]["labData"][i_nt]["sequenceData"]["files"]\
                              [index_read]["readOrder"] = key
            etl.wgs_submission_grz["donors"][0]["labData"][i_nt]["sequenceData"]["files"]\
                              [index_read]["readLength"] = gv.wgs_readLength

# add vcf info to submission_grz
with open(args.vcf_sha256_json, "r") as vcf_sha256:
    vcf_checksum = js.load(vcf_sha256)

with open(args.vcf_bytes_json, "r") as vcf_bytes:
   vcf_bytesize = js.load(vcf_bytes)

if vcf_checksum["vcf_checksum"][0]["file"] == vcf_bytesize["vcf_bytesize"][0]["file"]:
    vcf_name = vcf_checksum["vcf_checksum"][0]["file"]

vcf_sha256sum = vcf_checksum["vcf_checksum"][0]["fileChecksum"]

vcf_size = vcf_bytesize["vcf_bytesize"][0]["fileByteSize"]

index_tumor = 1
index_vcf = 2
if len(etl.wgs_submission_grz["donors"][0]["labData"][index_tumor]["sequenceData"]["files"]) == index_num_files[index_tumor]:
    # index 2 for vcf
    etl.wgs_submission_grz["donors"][0]["labData"][index_tumor]["sequenceData"]["files"]\
                      [index_vcf]["filePath"] = vcf_name
    etl.wgs_submission_grz["donors"][0]["labData"][index_tumor]["sequenceData"]["files"]\
                      [index_vcf]["fileType"] = gv.wgs_files_vcf
    etl.wgs_submission_grz["donors"][0]["labData"][index_tumor]["sequenceData"]["files"]\
                      [index_vcf]["checksumType"] = gv.wgs_checksumType
    etl.wgs_submission_grz["donors"][0]["labData"][index_tumor]["sequenceData"]["files"]\
                      [index_vcf]["fileChecksum"] = vcf_sha256sum
    etl.wgs_submission_grz["donors"][0]["labData"][index_tumor]["sequenceData"]["files"]\
                      [index_vcf]["fileSizeInBytes"] =  int(vcf_size)

#add bed info to submission_grz
#with open(args.bed_256_json, "r") as bed_256:
#    bed_checksum = js.load(bed_256)

#with open(args.bed_bytes_json, "r") as bed_bytes:
#    bed_bytesize = js.load(bed_bytes)

#if bed_checksum["bedfile_checksum"][0]["file"] == bed_bytesize["bed_bytesize"]\
#                                                           [0]["file"]:
#    bed_name = bed_checksum["bedfile_checksum"][0]["file"]

#bed_sha256sum = bed_checksum["bedfile_checksum"][0]["fileChecksum"]

#bed_size = bed_bytesize["bed_bytesize"][0]["fileByteSize"]

#index_bed = 3
#if len(etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]) == 4:

#    etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
#                      [index_bed]["filePath"] = rel_path_grz_cli+os.path.basename(bed_name)
#    etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
#                      [index_bed]["fileType"] = "bed"
#    etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
#                     [index_bed]["checksumType"] = "sha256"
#    etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
#                      [index_bed]["fileChecksum"] = bed_sha256sum
#    etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
#                      [index_bed]["fileSizeInBytes"] =  bed_size

# add cellularity
etl.wgs_submission_grz["donors"][0]["labData"][index_tumor]\
                  ["tumorCellCount"][0]\
                  ["count"] = int(wgs_final_page_dict["cellularity"])

etl.wgs_submission_grz["donors"][0]["labData"][index_tumor]\
                  ["tumorCellCount"][0]\
                  ["method"] = gv.wgs_tumorCellCount_method

#if patient_id == patient_id_check:

wgs_json = {
       "analysis" : p_data["analysis"],
       "mvtracker_uuid" : "",
       "pathoProId": p_data["pathoProId"],
       "nexusId": p_data["nexusId"],
       "molpathId": p_data["molpathId"],
       "patient_id" :  p_data["orbisId"],
       "firstname":p_data["firstname"],
       "lastame": p_data["lastname"],
       "dateofbrith": p_data["dateofbirth"],
       "network" : wgs_final_page_dict["network"],
       "diseaseType" : "oncological",
       "entity_excel" : wgs_final_page_dict["entity"],
       "entity_fm": p_data["entity"],
       "molecularReportod" : Oncology_Molecular_Report,
       "submission_grz" : etl.wgs_submission_grz
        }

   
with open(args.sample_id + "_submit.json", "w") as file:
    js.dump(
    wgs_json, 
    file, 
    indent=4, 
    sort_keys=False
)
        
