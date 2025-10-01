#!/usr/bin/env python

import pandas as pd
import os.path
import argparse
import json as js
import functions_pan_etl.MVdataset_generator_utils as etl
import functions_pan_etl.global_variables as gv
import requests
#import datetime

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
parser.add_argument("-p", "--patient_data_json", type=str)
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
    "tmb" : pan_final_page_dict["tmb"]["tmb_status"],
    "tmb_mutations_Mb" : pan_final_page_dict["tmb"]["tmb_mutations_Mb"],
    "msi_status" : pan_final_page_dict["msi"]["msi_status"],
    "msiAnzahlMarkerStabil" : pan_final_page_dict["msi"]["msiAnzahlMarkerStabil"],
    "msiAnzahlMarkerInstabil" : pan_final_page_dict["msi"]["msiAnzahlMarkerInstabil"],
    "msiAnzahlMarkerfail" : pan_final_page_dict["msi"]["msiAnzahlMarkerfail"],
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
etl.submission_grz["submission"]["submissionDate"] = "generated upon upload" #datetime.datetime.now().replace(hour=0, minute=0, second=0, microsecond=0).isoformat().split("T")[0]
etl.submission_grz["submission"]["submissionType"] = gv.submissionType # test modus
etl.submission_grz["submission"]["tanG"] = "data_from_MVtracker or MwTek"
etl.submission_grz["submission"]["localCaseId"] = "UUID???"
etl.submission_grz["submission"]["coverageType"] = "data_from_MVtracker or FileMaker"
etl.submission_grz["submission"]["submitterId"] = gv.submitterId
etl.submission_grz["submission"]["genomicDataCenterId"] = gv.genomicDataCenterId
etl.submission_grz["submission"]["clinicalDataNodeId"] = gv.clinicalDataNodeId
etl.submission_grz["submission"]["diseaseType"] = gv.diseaseType # assumed
etl.submission_grz["submission"]["genomicStudyType"] = gv.genomicStudyType # assumed
etl.submission_grz["submission"]["genomicStudySubtype"] = gv.genomicStudySubtype # assumed due to pancancer analysis
etl.submission_grz["submission"]["labName"] = gv.labName

# add broad consent information
etl.submission_grz["donors"][0]["donorPseudonym"] = gv.donorPseudonym
etl.submission_grz["donors"][0]["gender"] = p_data["gender"]
etl.submission_grz["donors"][0]["relation"] = gv.relation
etl.submission_grz["donors"][0]["mvConsent"]["presentationDate"] = "data_from_MVtracker"
etl.submission_grz["donors"][0]["mvConsent"]["version"] = "data_from_MVtracker"
etl.submission_grz["donors"][0]["mvConsent"]["scope"][0]["type"] = "data_from_MVtracker"
etl.submission_grz["donors"][0]["mvConsent"]["scope"][0]["date"] = "data_from_MVtracker"
etl.submission_grz["donors"][0]["mvConsent"]["scope"][0]["domain"] = "data_from_MVtracker"
etl.submission_grz["donors"][0]["researchConsents"][0]["schemaVersion"] = "data_from_MVtracker"
etl.submission_grz["donors"][0]["researchConsents"][0]["presentationDate"] = "data_from_MVtracker"
etl.submission_grz["donors"][0]["researchConsents"][0]["scope"] = "data_from_MVtracker"

# add info to sumbmission_grz - labData
etl.submission_grz["donors"][0]["labData"][0]\
                  ["labDataName"] = p_data["pathoProId"]

# add info tissue ontology
etl.submission_grz["donors"][0]["labData"][0]\
                  ["tissueOntology"]["name"] = gv.tissueOntology_name

etl.submission_grz["donors"][0]["labData"][0]\
                  ["tissueOntology"]["version"] = gv.tissueOntology_version

etl.submission_grz["donors"][0]["labData"][0]\
                  ["tissueTypeId"] = "data_from_FileMaker_Nexus_view"

etl.submission_grz["donors"][0]["labData"][0]\
                  ["tissueTypeName"] = "from icd03 xml bfarm or FileMaker_nexus_view"

# add further labData information

etl.submission_grz["donors"][0]["labData"][0]\
                  ["sampleDate"] = "date of the probe: FileMaker_nexus_view"

etl.submission_grz["donors"][0]["labData"][0]\
                  ["sampleConservation"] = "data_from_finalsheet_excel"

etl.submission_grz["donors"][0]["labData"][0]\
                  ["sequenceType"] = gv.pan_sequenceType # assumed due to pancancer analysis

etl.submission_grz["donors"][0]["labData"][0]\
                  ["sequenceSubtype"] = gv.pan_sequenceSubtype # assumed due to pancancer analysis

etl.submission_grz["donors"][0]["labData"][0]\
                  ["fragmentationMethod"] = gv.pan_fragmentationMethod #known

etl.submission_grz["donors"][0]["labData"][0]\
                  ["libraryType"] = "panel" # get from excel as well

etl.submission_grz["donors"][0]["labData"][0]\
                  ["libraryPrepKit"] = gv.pan_libraryPrepKit

etl.submission_grz["donors"][0]["labData"][0]\
                  ["libraryPrepKitManufacturer"] = gv.pan_libraryPrepKitManufacturer

etl.submission_grz["donors"][0]["labData"][0]\
                  ["sequencerModel"] = gv.pan_sequencerModel # or NextSeq550?

etl.submission_grz["donors"][0]["labData"][0]\
                  ["sequencerManufacturer"] = gv.pan_sequencerManufacturer # known

etl.submission_grz["donors"][0]["labData"][0]\
                  ["kitName"] = gv.pan_kitName

etl.submission_grz["donors"][0]["labData"][0]\
                  ["kitManufacturer"] = gv.pan_kitManufacturer

etl.submission_grz["donors"][0]["labData"][0]\
                  ["enrichmentKitManufacturer"] = gv.pan_enrichmentKitManufacturer

etl.submission_grz["donors"][0]["labData"][0]\
                  ["enrichmentKitDescription"] = gv.pan_enrichmentKitDescription

etl.submission_grz["donors"][0]["labData"][0]\
                  ["barcode"] = gv.pan_barcode

etl.submission_grz["donors"][0]["labData"][0]\
                  ["sequencingLayout"] = gv.pan_sequencingLayout # known

# add cellularity
etl.submission_grz["donors"][0]["labData"][0]\
                  ["tumorCellCount"][0]\
                  ["count"] = int(pan_final_page_dict["cellularity"])

etl.submission_grz["donors"][0]["labData"][0]\
                  ["tumorCellCount"][0]\
                  ["method"] = gv.pan_tumorCellCount_method

# add sequenceData info to grz_submission
etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]\
                  ["bioinformaticsPipelineName"] = gv.pan_bioinformaticsPipelineName

etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]\
                  ["bioinformaticsPipelineVersion"] = gv.pan_bioinformaticsPipelineVersion

etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]\
                  ["referenceGenome"] = gv.pan_referenceGenome

# add fastp qc info to submission_grz - "percentBasesAboveQualityThreshold"
with open(args.fastp_json, "r") as fastp_json:
    fastp_qc = js.load(fastp_json)

reads_q30_rate = round(fastp_qc["summary"]["before_filtering"]["q30_rate"],2) * 100

etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]\
                  ["percentBasesAboveQualityThreshold"] = round(float(reads_q30_rate),0)

# add bamfile info to submission_grz
with open(args.bam_json, "r") as bam_json:
    bam_qc = js.load(bam_json)

etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]\
                  ["meanDepthOfCoverage"] = round(float(bam_qc["bam_qc"][0]["mean_cov"]),2)

etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]\
                  ["minCoverage"] = round(float(bam_qc["bam_qc"][0]["min_cov"]),0)

etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]\
                  ["targetedRegionsAboveMinCoverage"] = float(bam_qc["bam_qc"][0]\
                                                        ["target_above_mincov"])

# add info of nonCodingVariants and caller information
etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]\
                  ["nonCodingVariants"] = gv.pan_nonCodingVariants

etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]\
                  ["callerUsed"][0]["name"] = gv.pan_callerUsed_name

etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]\
                  ["callerUsed"][0]["version"] = gv.pan_callerUsed_version

# the relative path of the grz_cli tool is always here sample_id/files/file
rel_path_grz_cli = args.sample_id+"/files/"

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
if len(etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]) == 4:
    # index 0 for read1
    # index 1 for read2
    for index_read in indices_reads:
        print(index_read)
        key, value = list(files_dict_cal_concat.items())[index_read]

        etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                          [index_read]["filePath"] = rel_path_grz_cli+files_dict_cal_concat[key]["file"]
        etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                          [index_read]["fileType"] = "fastq"
        etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                          [index_read]["checksumType"] = "sha256"
        etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                          [index_read]["fileChecksum"] = files_dict_cal_concat[key]["fileChecksum"]
        etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                          [index_read]["fileSizeInBytes"] = files_dict_cal_concat[key]["fileSizeInBytes"]
        etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
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
if len(etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]) == 4:
    # index 2 for vcf
    etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                      [index_vcf]["filePath"] = rel_path_grz_cli+vcf_name
    etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                      [index_vcf]["fileType"] = "vcf"
    etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                      [index_vcf]["checksumType"] = "sha256"
    etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                      [index_vcf]["fileChecksum"] = vcf_sha256sum
    etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
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
if len(etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]) == 4:

    etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                      [index_bed]["filePath"] = rel_path_grz_cli+os.path.basename(bed_name)
    etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                      [index_bed]["fileType"] = "bed"
    etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                      [index_bed]["checksumType"] = "sha256"
    etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                      [index_bed]["fileChecksum"] = bed_sha256sum
    etl.submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                      [index_bed]["fileSizeInBytes"] =  bed_size


#if patient_id == patient_id_check:

pancancer_json = {
       "analysis" : p_data["analysis"],
       "pathoProId": p_data["pathoProId"],
       "nexusId": p_data["nexusId"],
       "molpathId": p_data["molpathId"],
       "orbisId" :  p_data["orbisId"],
       "firstname":p_data["firstname"],
       "lastame": p_data["lastname"],
       "dateofbrith": p_data["dateofbirth"],
       "network" : "mv",
       "diseaseType" : "oncological",
       "entity_excel" : pan_final_page_dict["entity"],
       "entity_fm": p_data["entity"],
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
        
