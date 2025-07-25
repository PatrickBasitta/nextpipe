# -*- coding: utf-8 -*-
"""
@author: Patrick Basitta
"""

import os
import pandas as pd
import hashlib

def pan_snvdata_to_dicts(filepath,excel_file,patient_id,idx):
    
    with open(filepath, 'rb') as snv_file: 
        snv_data = pd.read_excel(snv_file, 
                                 sheet_name= excel_file.sheet_names[idx],
                                 engine="openpyxl")   
        
        # get idx/IDX of reported x/X
        idx_of_variants_to_report = snv_data.index\
                                   [snv_data["submit"] == "x"].tolist()
                                   
        IDX_of_variants_to_report = snv_data.index\
                                   [snv_data["submit"] == "X"].tolist()
                                   
        # merge idx/IDX and sort
        merged = idx_of_variants_to_report + IDX_of_variants_to_report
        merged.sort()

        # get variants                                
        variants_to_report = snv_data.loc[merged, :]
        
        # get all necessary columns !!!change pancancer python script!!!
        
        # add column identifier, genomicSource = "somatic" and LOH
        variants_to_report["identifier"] = ""
        variants_to_report["genomicSource"] = "somatic"
        variants_to_report["LOH"] = "not_available"
        
        # rename colums
        variants_to_report = variants_to_report.rename(columns = 
                                       {"Gen":"gene.diplay",
                                        "HGNC_MV":"gene.code",
                                        "Transcript_ID": "transcriptId",
                                        "cDNA_Change": "dnaChange",
                                        "Amino_Acid_Change": "proteinChange",
                                        "Consequence": "variantTypes"
                                        })
        
        # get final columns to report
        final_data = variants_to_report[["identifier",
                                         "genomicSource",
                                         "gene.diplay",
                                         "gene.code",
                                         "transcriptId",
                                         "dnaChange",
                                         "proteinChange",
                                         "localization",
                                         "variantTypes",
                                         "LOH"]]
        
        # VariantID-generator
        smallVariantId_lst = []
        for i in range(len(final_data)):
            i = i + 1
            smallVariantId_lst.append("smallVariantId_"+str(i))
            
        # add id to identifier column
        final_data.loc[:,"identifier"] =  smallVariantId_lst
        
        # make dictionaries
        smallVariatns_records = final_data.to_dict("records")
        
        return smallVariatns_records
    

def pan_final_page_to_dict(filepath,excel_file,idx1):
    
        with open(filepath, 'rb') as file: 
            pancancer_page_final = pd.read_excel(file, 
                                     sheet_name= excel_file.sheet_names[idx1],
                                     engine="openpyxl")

            report_dict = dict()
            
            # entity
            # entity_dict = dict()
            report_dict["entity"] = pancancer_page_final.loc[0, "Entitaet"]   
            
            # cellularity
            report_dict["cellularity"] = pancancer_page_final.loc[0, "TZ"] 
            
            # Average_coverage
            report_dict["average_coverage"] = pancancer_page_final.loc[0, "Average_coverage"] 
            
            # Coverage>100
            report_dict["overage_above100"] = pancancer_page_final.loc[0, "Coverage>100"]
            
            # MSI as dict
            msi = dict()
            msi["msi_status"] = pancancer_page_final.loc[0, "MSI"]
            msi["msiAnzahlMarkerStabil"] = int(pancancer_page_final.loc\
                                          [0, "Anzahl_Marker_stabil"])
            msi["msiAnzahlMarkerInstabil"] = int(pancancer_page_final.loc\
                                          [0, "Anzahl_Marker_instabil"])
            msi["msiAnzahlMarkerfail"] = int(pancancer_page_final.loc\
                                          [0, "Anzahl_Marker_fail"])
            report_dict["msi"] = msi    
            # TMB
            tmb = dict()
            tmb["tmb_status"] =  pancancer_page_final.loc[0, "TMB"]
            tmb["tmb_mutations_Mb"] =  int(pancancer_page_final.loc[0, "mutations_Mb"])
            report_dict["tmb"] = tmb 
            
            return report_dict
                    

submission_grz = {
    "submissionDate" : "",
    "submissionType" : "",
    "tanG": "",
    "localCaseId" : "",
    "coverageType" : "",
    "submitterId" : "",
    "genomicDataCenterId" : "",
    "clinicalDataNodeId" : "",
    "diseaseType" : "",
    "genomicStudyType" : "",
    "genomicStudySubtype" : "",
    "labName" : "",
    "donors" : [{
        "donorPseudonym" : "",
        "gender" : "",
        "relation" : ""
        }],
    "mvConsent" : {
        "presentationDate" : "",
        "version" : ""
        },
    "scope" : [{
        "type" : "",
        "date" : "",
        "domain" : ""
        }],
    "researchConsents" : [{
        "schemaVersion" : "",
        "presentationDate" : "",
        "scope" : ""
        }],
    "labData" : [{
        "labDataName" : "FFPE DNA tumor",
        "tissueOntology" : {
            "name" : "not_yet_available",
            "version" : "not_yet_available"
            },
        "tissueTypeId" : "not_yet_available",
        "tissueTypeName" : "not_yet_available",
        "sampleDate" : "", 
        "sampleConservation" : "",
        "sequenceType" : "", 
        "sequenceSubtype" : "not_yet_available",
        "fragmentationMethod" : "not_yet_available",
        "libraryType" : "panel",
        "libraryPrepKit" : "not_yet_available",
        "libraryPrepKitManufacturer" : "not_yet_available",
        "sequencerModel" : "not_yet_available",
        "sequencerManufacturer" : "not_yet_available",
        "kitName" : "not_yet_available",
        "kitManufacturer" : "not_yet_available",
        "enrichmentKitManufacturer" : "not_yet_available",
        "enrichmentKitDescription" : "not_yet_available",
        "barcode" : "na",
        "sequencingLayout" : "paired-end",
        "tumorCellCount" : [{
            "count" : "not_yet_available",
            "method" : "pathology"
            }],       
        }],
    "sequenceData" : {
        "bioinformaticsPipelineName" : "not_yet_available", 
        "bioinformaticsPipelineVersion" : "not_yet_available",
        "referenceGenome" : "GRCh38",
        "percentBasesAboveQualityThreshold" : "not_yet_available",
        "meanDepthOfCoverage" : "not_yet_available",
        "minCoverage" : "not_yet_available",
        "targetedRegionsAboveMinCoverage" : "not_yet_available",
        "files" : "not_yet_available",
        "nonCodingVariants" : "not_yet_available",
        "callerUsed" : [{
            "name" : "not_yet_available",
            "version" : "not_yet_available"
            }],
        },
    "files" : [{
        "filePath" : "not_yet_available",
        "fileType" : "fastq",
        "checksumType" :"sha256",
        "fileChecksum" : "not_yet_available",
        "fileSizeInBytes" : "not_yet_available",
        "readOrder" : "R1",
        "readLenght" : "not_yet_available",
        "flowcellId" : "not_yet_available",
        "laneId" : "not_yet_available"
        },
        {
        "filePath" : "not_yet_available",
        "fileType" : "fastq",
        "checksumType" :"sha256",
        "fileChecksum" : "not_yet_available",
        "fileSizeInBytes" : "not_yet_available",
        "readOrder" : "R2",
        "readLenght" : "not_yet_available",
        "flowcellId" : "not_yet_available",
        "laneId" : "not_yet_available"
        },
        {
        "filePath" : "not_yet_available",
        "fileType" : "vcf",
        "checksumType" :"sha256",
        "fileChecksum" : "not_yet_available",
        "fileSizeInBytes" : "not_yet_available"
        },
        {
        "filePath" : "not_yet_available",
        "fileType" : "not_yet_available",
        "checksumType" :"sha256",
        "fileChecksum" : "not_yet_available",
        "fileSizeInBytes" : "not_yet_available"
        }]
    }

#do files calculations
#def sha256sum(fastq_file):
#    with open(fastq_file, "rb", buffering=0) as fq_file:
#        return hashlib.file_digest(fq_file, "sha256").hexdigest()

def get_fastq_info(lst_R1_R2_dict,number_reads):
    collect_fq_file_name = []
    #if os.path.exists(grz_cli_dir):
    for index_read_fq in range(number_reads):
        for key in lst_R1_R2_dict[index_read_fq].keys():
            print(key)
                
            for fastq_file in lst_R1_R2_dict[index_read_fq][key].values():
                print(fastq_file)
                
            fq_file_name = os.path.basename(fastq_file)
            collect_fq_file_name.append(fq_file_name)
    
            # calulate sha256sum and store in dict
            #os.path.basename(fastq_file) use grz_cli path!!! or compare sums!!!
            #with open(fastq_file, "rb", buffering=0) as fq_file:
            #    data = fq_file.read()
            #    hashlib.update(data)
            #    digest = hashlib.hexdigest()
                #sha256sum = hashlib.file_digest(fq_file, "sha256").hexdigest()
            #sha256sum_fq_file = digest
            #add_sha256sum = {"fileChecksum": sha256sum_fq_file}
            #lst_R1_R2_dict[index_read_fq][key].update(add_sha256sum)
          
            # get file size in bytes
            file_size_in_bytes = os.path.getsize(fastq_file)
            add_fileSizeBytes = {"fileSizeInBytes": file_size_in_bytes}
            lst_R1_R2_dict[index_read_fq][key].update(add_fileSizeBytes)
                
    return collect_fq_file_name, lst_R1_R2_dict
