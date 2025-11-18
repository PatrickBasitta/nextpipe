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
"SmallVariants" : smallVariatns_records,
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
etl.pan_submission_grz["submission"]["submissionDate"] = "generated upon upload" #datetime.datetime.now().replace(hour=0, minute=0, second=0, microsecond=0).isoformat().split("T")[0]
etl.pan_submission_grz["submission"]["submissionType"] = gv.submissionType
etl.pan_submission_grz["submission"]["tanG"] = ""
etl.pan_submission_grz["submission"]["localCaseId"] = ""
etl.pan_submission_grz["submission"]["coverageType"] = gv.coverageType
etl.pan_submission_grz["submission"]["submitterId"] = gv.submitterId
etl.pan_submission_grz["submission"]["genomicDataCenterId"] = gv.genomicDataCenterId
etl.pan_submission_grz["submission"]["clinicalDataNodeId"] = gv.clinicalDataNodeId
etl.pan_submission_grz["submission"]["diseaseType"] = gv.diseaseType # assumed
etl.pan_submission_grz["submission"]["genomicStudyType"] = gv.genomicStudyType # assumed
etl.pan_submission_grz["submission"]["genomicStudySubtype"] = gv.pan_genomicStudySubtype # assumed due to pancancer analysis
etl.pan_submission_grz["submission"]["labName"] = gv.labName

# add broad consent information
etl.pan_submission_grz["donors"][0]["donorPseudonym"] = gv.pan_donorPseudonym
etl.pan_submission_grz["donors"][0]["gender"] = p_data["gender"]
etl.pan_submission_grz["donors"][0]["relation"] = gv.pan_relation
#etl.pan_submission_grz["donors"][0]["mvConsent"]["presentationDate"] = ""
#etl.pan_submission_grz["donors"][0]["mvConsent"]["version"] = ""
#etl.pan_submission_grz["donors"][0]["mvConsent"]["scope"][0]["type"] = ""
#etl.pan_submission_grz["donors"][0]["mvConsent"]["scope"][0]["date"] = ""
#etl.pan_submission_grz["donors"][0]["mvConsent"]["scope"][0]["domain"] = ""
#etl.pan_submission_grz["donors"][0]["researchConsents"][0]["schemaVersion"] = ""
#etl.pan_submission_grz["donors"][0]["researchConsents"][0]["presentationDate"] = ""
#etl.pan_submission_grz["donors"][0]["researchConsents"][0]["scope"] = ""


# add info to sumbmission_grz - labData
etl.pan_submission_grz["donors"][0]["labData"][0]\
                  ["labDataName"] = gv.pan_labDataName

# add info tissue ontology
etl.pan_submission_grz["donors"][0]["labData"][0]\
                  ["tissueOntology"]["name"] = gv.pan_tissueOntology_name

etl.pan_submission_grz["donors"][0]["labData"][0]\
                  ["tissueOntology"]["version"] = gv.pan_tissueOntology_version

etl.pan_submission_grz["donors"][0]["labData"][0]\
                  ["tissueTypeId"] = ""

etl.pan_submission_grz["donors"][0]["labData"][0]\
                  ["tissueTypeName"] = ""

# add further labData information

etl.pan_submission_grz["donors"][0]["labData"][0]\
                  ["sampleDate"] = pan_final_page_dict["sampledate_T"]

etl.pan_submission_grz["donors"][0]["labData"][0]\
                  ["sampleConservation"] = pan_final_page_dict["sampleConservation_T"]

etl.pan_submission_grz["donors"][0]["labData"][0]\
                  ["sequenceType"] = gv.pan_sequenceType # assumed due to pancancer analysis

etl.pan_submission_grz["donors"][0]["labData"][0]\
                  ["sequenceSubtype"] = gv.pan_sequenceSubtype # assumed due to pancancer analysis

etl.pan_submission_grz["donors"][0]["labData"][0]\
                  ["fragmentationMethod"] = gv.pan_fragmentationMethod #known

etl.pan_submission_grz["donors"][0]["labData"][0]\
                  ["libraryType"] = gv.pan_libraryType # get from excel as well

etl.pan_submission_grz["donors"][0]["labData"][0]\
                  ["libraryPrepKit"] = gv.pan_libraryPrepKit

etl.pan_submission_grz["donors"][0]["labData"][0]\
                  ["libraryPrepKitManufacturer"] = gv.pan_libraryPrepKitManufacturer

etl.pan_submission_grz["donors"][0]["labData"][0]\
                  ["sequencerModel"] = pan_final_page_dict["sequencer"] # or NextSeq550?

etl.pan_submission_grz["donors"][0]["labData"][0]\
                  ["sequencerManufacturer"] = gv.pan_sequencerManufacturer # known

etl.pan_submission_grz["donors"][0]["labData"][0]\
                  ["kitName"] = pan_final_page_dict["kit_name"]

etl.pan_submission_grz["donors"][0]["labData"][0]\
                  ["kitManufacturer"] = gv.pan_kitManufacturer

etl.pan_submission_grz["donors"][0]["labData"][0]\
                  ["enrichmentKitManufacturer"] = gv.pan_enrichmentKitManufacturer

etl.pan_submission_grz["donors"][0]["labData"][0]\
                  ["enrichmentKitDescription"] = gv.pan_enrichmentKitDescription

etl.pan_submission_grz["donors"][0]["labData"][0]\
                  ["barcode"] = pan_final_page_dict["barcode_T"]

etl.pan_submission_grz["donors"][0]["labData"][0]\
                  ["sequencingLayout"] = gv.pan_sequencingLayout # known

# add cellularity
etl.pan_submission_grz["donors"][0]["labData"][0]\
                  ["tumorCellCount"][0]\
                  ["count"] = int(pan_final_page_dict["cellularity"])

etl.pan_submission_grz["donors"][0]["labData"][0]\
                  ["tumorCellCount"][0]\
                  ["method"] = gv.pan_tumorCellCount_method

# add sequenceData info to grz_submission
etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]\
                  ["bioinformaticsPipelineName"] = gv.pan_bioinformaticsPipelineName

etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]\
                  ["bioinformaticsPipelineVersion"] = gv.pan_bioinformaticsPipelineVersion

etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]\
                  ["referenceGenome"] = gv.pan_referenceGenome

# add fastp qc info to submission_grz - "percentBasesAboveQualityThreshold"
with open(args.fastp_json, "r") as fastp_json:
    fastp_qc = js.load(fastp_json)

etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]\
                      ["percentBasesAboveQualityThreshold"]["minimumQuality"] = gv.pan_minimumQuality

reads_q30_rate = round(fastp_qc["summary"]["before_filtering"]["q30_rate"],2) * 100

etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]\
                  ["percentBasesAboveQualityThreshold"]["percent"] = round(float(reads_q30_rate),0)

# add bamfile info to submission_grz
with open(args.bam_json, "r") as bam_json:
    bam_qc = js.load(bam_json)

etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]\
                  ["meanDepthOfCoverage"] = round(float(bam_qc["bam_qc"][0]["mean_cov"]),2)

etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]\
                  ["minCoverage"] = gv.pan_minCoverage

etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]\
                  ["targetedRegionsAboveMinCoverage"] = float(bam_qc["bam_qc"][0]\
                                                        ["targets_above_mincov"])

# add info of nonCodingVariants and caller information
etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]\
                  ["nonCodingVariants"] = gv.pan_nonCodingVariants

etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]\
                  ["callerUsed"][0]["name"] = gv.pan_callerUsed_name

etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]\
                  ["callerUsed"][0]["version"] = gv.pan_callerUsed_version

# the relative path of the grz_cli tool is always here sample_id/files/file
#rel_path_grz_cli = args.sample_id+"/files/"

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
if len(etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]) == 4:
    # index 0 for read1
    # index 1 for read2
    for index_read in indices_reads:
        print(index_read)
        key, value = list(files_dict_cal_concat.items())[index_read]

        etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                          [index_read]["filePath"] = files_dict_cal_concat[key]["file"]
        etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                          [index_read]["fileType"] = gv.pan_files_fastq
        etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                          [index_read]["checksumType"] = gv.pan_checksumType
        etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                          [index_read]["fileChecksum"] = files_dict_cal_concat[key]["fileChecksum"]
        etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                          [index_read]["fileSizeInBytes"] = int(files_dict_cal_concat[key]["fileSizeInBytes"])
        etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                          [index_read]["readOrder"] = key
        etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                          [index_read]["readLength"] = gv.pan_readLength


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
if len(etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]) == 4:
    # index 2 for vcf
    etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                      [index_vcf]["filePath"] = vcf_name
    etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                      [index_vcf]["fileType"] = gv.pan_files_vcf
    etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                      [index_vcf]["checksumType"] = gv.pan_checksumType
    etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                      [index_vcf]["fileChecksum"] = vcf_sha256sum
    etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                      [index_vcf]["fileSizeInBytes"] =  int(vcf_size)

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
if len(etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]) == 4:

    etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                      [index_bed]["filePath"] = os.path.basename(bed_name)
    etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                      [index_bed]["fileType"] = gv.pan_files_bed
    etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                      [index_bed]["checksumType"] = gv.pan_checksumType
    etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                      [index_bed]["fileChecksum"] = bed_sha256sum
    etl.pan_submission_grz["donors"][0]["labData"][0]["sequenceData"]["files"]\
                      [index_bed]["fileSizeInBytes"] =  int(bed_size)


#if patient_id == patient_id_check:

pancancer_json = {
       "analysis" : p_data["analysis"],
       "mvtracker_uuid" : "",
       "pathoProId": p_data["pathoProId"],
       "nexusId": p_data["nexusId"],
       "molpathId": p_data["molpathId"],
       "patient_id" :  p_data["orbisId"],
       "firstname":p_data["firstname"],
       "lastame": p_data["lastname"],
       "dateofbrith": p_data["dateofbirth"],
       "network" : pan_final_page_dict["network"],
       "diseaseType" : "oncological",
       "entity_excel" : pan_final_page_dict["entity"],
       "entity_fm": p_data["entity"],
       "molecularReportod" : Oncology_Molecular_Report,
       "submission" : etl.pan_submission_grz
        }
    
with open(args.sample_id + "_submit.json", "w") as file:
    js.dump(
    pancancer_json, 
    file, 
    indent=4, 
    sort_keys=False
)
        
