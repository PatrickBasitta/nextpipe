# -*- coding: utf-8 -*-
"""
@author: Patrick Basitta
"""

import os
import pandas as pd
import hashlib
import requests
import json as js
import functions_etl.global_variables as gv
import xml.etree.ElementTree as ET

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
        #variants_to_report["interpretation"] = ""
        variants_to_report["LOH"] = "not_available"
        
        # rename colums
        variants_to_report = variants_to_report.rename(columns = 
                                       {"Chromosom":"chromosome",
                                        "Position":"startPosition",
                                        "End Position": "endPosition",
                                        "Reference": "ref",
                                        "Allele": "alt",
                                        "Gen":"gene.diplay",
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
                                         "chromosome",
                                         "startPosition",
                                         "endPosition",
                                         "ref",
                                         "alt",
                                         "localization",
                                         "variantTypes",
                                         "interpretation",
                                         "LOH"]]
        
        final_data["interpretation"] = final_data["interpretation"].apply(lambda x: int(x))

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

def wxs_snvdata_to_dicts(filepath,excel_file,patient_id,idx,hgnc):
    
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
        
        # get all necessary columns
        
        # add column identifier, genomicSource = "somatic" and LOH
        variants_to_report["identifier"] = ""
        variants_to_report["genomicSource"] = "somatic"
        #variants_to_report["interpretation"] = ""
        variants_to_report["LOH"] = "not_available"
        
        # rename colums
        variants_to_report = variants_to_report.rename(columns = 
                                       {"Chromosome":"chromosome",
                                        "Start_Position":"startPosition",
                                        "End_Position": "endPosition",
                                        "Reference_Allele": "ref",
                                        "Tumor_Seq_Allele2": "alt",
                                        "HUGO_SYMBOL":"gene.display",
                                        "NM-Nummer": "transcriptId",
                                        "HGVSc": "dnaChange",
                                        "HGVSp": "proteinChange",
                                        "CSQ_Consequence": "variantTypes"
                                        })
        
        # add localization, due not present in xlsx (in PanCancer it is present)
        
        # getting fields {coding,inRegulatoryElements,
        # notInCodingAndNotInRegulatoryElements} in field localization
        # required for the MV Oncology report
        variants_to_report["localization"] = "-"

        # field: coding (bool: according to ensemble variation catalog and HGVSp = xxx.p.xxx / exists)
        # field: notInCodingAndNotInRegulatoryElements (bool: )
        # field: inRegulatoryElements(bool: )
        
        # Reset indices
        variants_to_report = variants_to_report.reset_index(drop="TRUE")
        
        # coding_SO_term
        coding_SO_term = ["synonymous_variant",
                          "missense_variant", 
                          "inframe_insertion", 
                          "inframe_deletion", 
                          "stop_gained",
                          "frameshift_variant", 
                          "coding_sequence_variant",
                          "start_lost", 
                          "stop_lost", 
                          "start_retained_variant",
                          "stop_retained_variant", 
                          "incomplete_terminal_codon_variant",
                          "transcript_ablation",
                          "transcript_amplification", 
                          "protein_altering_variant",
                          "coding_transcript_variant",
                          "NMD_transcript_variant"]
        
        # regulatroy_SO_term
        regulatroy_SO_term = ["TFBS_ablation", 
                              "TFBS_amplification",
                              "TF_binding_site_variant", 
                              "regulatory_region_ablation",
                              "regulatory_region_amplification",
                              "regulatory_region_variant"]
        
        # define "localization"
        for i in range(len(variants_to_report["proteinChange"])):
            if variants_to_report.loc[i,"variantTypes"] in coding_SO_term and pd.isna\
               (variants_to_report.loc[i,"proteinChange"]) == False and variants_to_report.loc\
                     [i,"proteinChange"].startswith("p."):
                         variants_to_report.loc[i,"localization"] = "coding"            
            elif variants_to_report.loc\
                [i,"variantTypes"] in regulatroy_SO_term and variants_to_report.loc\
                      [i,"proteinChange"].startswith("p.") == False:
                          variants_to_report.loc[i,"localization"] = "inRegulatoryElements"
            elif variants_to_report.loc\
                [i,"variantTypes"] not in regulatroy_SO_term and variants_to_report.loc\
                    [i,"variantTypes"] not in coding_SO_term and str(variants_to_report.loc\
                      [i,"proteinChange"]).startswith("p.") == False:
                          variants_to_report.loc[i,"localization"] = "notInCodingAndNotInRegulatoryElements"

        # fields {inRegulatoryElements,notInCodingAndNotInRegulatoryElements} are
        # preliminary not catched here in PanCancer, but must be added in WES/WGS
        
        # add HGNC_MV, due not present in xlsx (in PanCancer it is present)
        # make HGNC_code column
        hgnc["HGNC_code"] = hgnc["HGNC"]
        # make dict (sorts automatically)
        hgnc_dict = dict(zip(hgnc["ID"] ,hgnc["HGNC_code"]))
        # add hgnc code
        for i in range(len(variants_to_report["gene.display"])):
            variants_to_report.loc[i, "gene.code"] = hgnc_dict[variants_to_report.loc[i, "gene.display"]]
        
        # convert np int64 to int
        variants_to_report["interpretation_int"] = ""
        for number in range(len(variants_to_report["interpretation"])):
            variants_to_report.loc[number,"interpretation_int"] = int(variants_to_report.loc\
                                                                     [number,"interpretation"])
        # get final columns to report
        final_data = variants_to_report[["identifier",
                                         "genomicSource",
                                         "gene.display",
                                         "gene.code",
                                         "transcriptId",
                                         "dnaChange",
                                         "proteinChange",
                                         "chromosome",
                                         "startPosition",
                                         "endPosition",
                                         "ref",
                                         "alt",
                                         "localization",
                                         "variantTypes",
                                         "interpretation_int",
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
            report_dict["coverage_above100"] = pancancer_page_final.loc[0, "Coverage>100"]
            
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

def wxs_final_page_to_dict(filepath,excel_file,idx1):
    
        with open(filepath, 'rb') as file: 
            wxs_page_final = pd.read_excel(file, 
                                     sheet_name= excel_file.sheet_names[idx1],
                                     engine="openpyxl")

            report_dict = dict()
            
            # network
            report_dict["network"] = wxs_page_final.loc[0,"network"]
            
            # library_type
            report_dict["library_type"] = wxs_page_final.loc[0,"library_type"]
            
            # sampleConservation_T
            report_dict["sampleConservation_T"] = wxs_page_final.loc[0,"sampleConservation_T"]
            
            # sampleConservation_N
            report_dict["sampleConservation_N"] = wxs_page_final.loc[0,"sampleConservation_N"]
                        
            # entity
            # entity_dict = dict()
            # report_dict["entity"] = wxs_page_final.loc[0, "Entitaet"]   
            
            # sampledates
            report_dict["sampledate_T"] = wxs_page_final.loc[0, "sampledate_T"]

            report_dict["sampledate_N"] = wxs_page_final.loc[0, "sampledate_N"]

            # kit_name
            report_dict["kit_name"] = wxs_page_final.loc[0, "kit_name"]

            # barcodes
            report_dict["barcode_T"] = wxs_page_final.loc[0, "barcode_T"]

            report_dict["barcode_N"] = wxs_page_final.loc[0, "barcode_N"]

            # sequencer
            report_dict["sequencer"] = wxs_page_final.loc[0, "sequencer"]

            # entity
            # entity_dict = dict()
            report_dict["entity"] = wxs_page_final.loc[0, "Entitaet"]

            # cellularity
            report_dict["cellularity"] = float(wxs_page_final.loc[0, "TZ (%)"]) 
            
            # Average_coverage !
            #report_dict["average_coverage"] = wxs_page_final.loc[0, "Average_coverage"] 
            
            # Coverage>100 !
            #report_dict["overage_above100"] = wxs_page_final.loc[0, "Coverage>100"]
            
            # MSI as dict
            msi = dict()
            msi["msi_status"] = wxs_page_final.loc[0, "MSI_UKB"]
            msi["Ergebnis_MSIsensor_pro"] = float(wxs_page_final.loc\
                                          [0, "Ergebnis_MSIsensor_pro(%)_UKB"])
            msi["msi_OA"] = wxs_page_final.loc[0, "MSI_OA"]

            report_dict["msi"] = msi   
            
            # TMB
            tmb = dict()
            tmb["tmb_status"] =  wxs_page_final.loc[0, "TMB_UKB"]
            tmb["Anzahl_Mutationen_missense"] =  wxs_page_final.loc\
                                                 [0, "Anzahl_Mutationen_missense_UKB"]
            tmb["mutations_Mb"] =  wxs_page_final.loc\
                                                 [0, "mutations_Mb_UKB"]
            tmb["tmb_OA"] = wxs_page_final.loc[0, "TMB_OA"]
            tmb["mutations_Mb_OA"] = wxs_page_final.loc[0, "mutations_Mb_OA"]
             
            report_dict["tmb"] = tmb

            # HRD
            hrd = dict()
            hrd = wxs_page_final.loc[0, "HR deficiency score_OA"]
            
            return report_dict
                    

pan_submission_grz = {
    "submission" : {
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
        "labName" : ""
         },
    "donors" : [{
        "donorPseudonym" : "",
        "gender" : "",
        "relation" : "",
        "mvConsent" : {
            "presentationDate" : "",
            "version" : "",
            "scope" : [{
                "type" : "",
                "date" : "",
                "domain" : ""
                }],
            },
        "researchConsents" : [{
            "schemaVersion" : "",
            "presentationDate" : "",
            "scope" : ""
            }],
        "labData" : [{
            "labDataName" : "",
            "tissueOntology" : {
                "name" : "",
                "version" : ""
                },
            "tissueTypeId" : "",
            "tissueTypeName" : "",
            "sampleDate" : "",
            "sampleConservation" : "",
            "sequenceType" : "",
            "sequenceSubtype" : "",
            "fragmentationMethod" : "",
            "libraryType" : "",
            "libraryPrepKit" : "",
            "libraryPrepKitManufacturer" : "",
            "sequencerModel" : "",
            "sequencerManufacturer" : "",
            "kitName" : "",
            "kitManufacturer" : "",
            "enrichmentKitManufacturer" : "",
            "enrichmentKitDescription" : "",
            "barcode" : "",
            "sequencingLayout" : "",
            "tumorCellCount" : [{
                "count" : "",
                "method" : ""
                }],
            "sequenceData" : {
                "bioinformaticsPipelineName" : "",
                "bioinformaticsPipelineVersion" : "",
                "referenceGenome" : "",
                "percentBasesAboveQualityThreshold" : "",
                "meanDepthOfCoverage" : "",
                "minCoverage" : "",
                "targetedRegionsAboveMinCoverage" : "",
                "nonCodingVariants" : "",
                "callerUsed" : [{
                    "name" : "",
                    "version" : ""
                    }],
                "files" : [{
                    "filePath" : "",
                    "fileType" : "",
                    "checksumType" :"",
                    "fileChecksum" : "",
                    "fileSizeInBytes" : "",
                    "readOrder" : "",
                    #"readLenght" : "",
                    #"flowcellId" : "",
                    #"laneId" : ""
                    },
                    {
                    "filePath" : "",
                    "fileType" : "",
                    "checksumType" :"",
                    "fileChecksum" : "",
                    "fileSizeInBytes" : "",
                    "readOrder" : "",
                    #"readLenght" : "",
                    #"flowcellId" : "",
                    #"laneId" : ""
                    },
                    {
                    "filePath" : "",
                    "fileType" : "",
                    "checksumType" :"",
                    "fileChecksum" : "",
                    "fileSizeInBytes" : ""
                    },
                    {
                    "filePath" : "",
                    "fileType" : "",
                    "checksumType" :"",
                    "fileChecksum" : "",
                    "fileSizeInBytes" : ""
                    }],
               },
          }],
    }],
}

wes_submission_grz = {
    "$schema": "https://raw.githubusercontent.com/BfArM-MVH/MVGenomseq/refs/tags/v1.2.1/GRZ/grz-schema.json",
    "submission" : {
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
        "labName" : ""
         },
    "donors" : [{
        "donorPseudonym" : "",
        "gender" : "",
        "relation" : "",
        "mvConsent" : {
            #"presentationDate" : "",
            "version" : "",
            "scope" : [{
                "type" : "",
                "date" : "",
                "domain" : "mvSequencing"
                },
                {
                "type" : "",
                "date" : "",
                "domain" : "reIdentification"
                },
                {
                "type" : "",
                "date" : "",
                "domain" : "caseIdentification"
                }
            ],
        },
        "researchConsents" : [
        ],
        #"researchConsents" : [{
        #    "schemaVersion" : "",
        #    "presentationDate" : "",
        #    "scope" : ""
        #    }],
        "labData" : [{
            "labDataName" : "",
            "tissueOntology" : {
                "name" : "",
                "version" : ""
                },
            "tissueTypeId" : "",
            "tissueTypeName" : "",
            "sampleDate" : "",
            "sampleConservation" : "",
            "sequenceType" : "",
            "sequenceSubtype" : "",
            "fragmentationMethod" : "",
            "libraryType" : "",
            "libraryPrepKit" : "",
            "libraryPrepKitManufacturer" : "",
            "sequencerModel" : "",
            "sequencerManufacturer" : "",
            "kitName" : "",
            "kitManufacturer" : "",
            "enrichmentKitManufacturer" : "",
            "enrichmentKitDescription" : "",
            "barcode" : "",
            "sequencingLayout" : "",
            "sequenceData" : {
                "bioinformaticsPipelineName" : "",
                "bioinformaticsPipelineVersion" : "",
                "referenceGenome" : "",
                "percentBasesAboveQualityThreshold" : {
                    "minimumQuality": "",
                    "percent": ""
                 },
                "meanDepthOfCoverage" : "",
                "minCoverage" : "",
                "targetedRegionsAboveMinCoverage" : "",
                "nonCodingVariants" : "",
                "callerUsed" : [{
                    "name" : "",
                   "version" : ""
                    }],
                "files" : [{
                    "filePath" : "",
                    "fileType" : "",
                    "checksumType" :"",
                    "fileChecksum" : "",
                    "fileSizeInBytes" : "",
                    "readOrder" : "",
                    "readLength" : "",
                    #"flowcellId" : "",
                    #"laneId" : ""
                    },
                    {
                    "filePath" : "",
                    "fileType" : "",
                    "checksumType" :"",
                    "fileChecksum" : "",
                    "fileSizeInBytes" : "",
                    "readOrder" : "",
                    "readLength" : "",
                    #"flowcellId" : "",
                    #"laneId" : ""
                    },
                    {
                    "filePath" : "",
                    "fileType" : "",
                    "checksumType" :"",
                    "fileChecksum" : "",
                    "fileSizeInBytes" : ""
                    }],
               },
          },
          {
            "labDataName" : "",
            "tissueOntology" : {
                "name" : "",
                "version" : ""
                },
            "tissueTypeId" : "",
            "tissueTypeName" : "",
            "sampleDate" : "",
            "sampleConservation" : "",
            "sequenceType" : "",
            "sequenceSubtype" : "",
            "fragmentationMethod" : "",
            "libraryType" : "",
            "libraryPrepKit" : "",
            "libraryPrepKitManufacturer" : "",
            "sequencerModel" : "",
            "sequencerManufacturer" : "",
            "kitName" : "",
            "kitManufacturer" : "",
            "enrichmentKitManufacturer" : "",
            "enrichmentKitDescription" : "",
            "barcode" : "",
            "sequencingLayout" : "",
            "tumorCellCount" : [{
                "count" : "",
                "method" : ""
                }],
            "sequenceData" : {
                "bioinformaticsPipelineName" : "",
                "bioinformaticsPipelineVersion" : "",
                "referenceGenome" : "",
                "percentBasesAboveQualityThreshold" : {
                     "minimumQuality": "",
                     "percent": ""
                },
                "meanDepthOfCoverage" : "",
                "minCoverage" : "",
                "targetedRegionsAboveMinCoverage" : "",
                "nonCodingVariants" : "",
                "callerUsed" : [{
                    "name" : "",
                    "version" : ""
                    }],
                "files" : [{
                    "filePath" : "",
                    "fileType" : "",
                    "checksumType" :"",
                    "fileChecksum" : "",
                    "fileSizeInBytes" : "",
                    "readOrder" : "",
                    "readLength" : "",
                    #"flowcellId" : "",
                    #"laneId" : ""
                    },
                    {
                    "filePath" : "",
                    "fileType" : "",
                    "checksumType" :"",
                    "fileChecksum" : "",
                    "fileSizeInBytes" : "",
                    "readOrder" : "",
                    "readLength" : "",
                    #"flowcellId" : "",
                    #"laneId" : ""
                    },
                    {
                    "filePath" : "",
                    "fileType" : "",
                    "checksumType" :"",
                    "fileChecksum" : "",
                    "fileSizeInBytes" : ""
                    },
                     {
                    "filePath" : "",
                    "fileType" : "",
                    "checksumType" :"",
                    "fileChecksum" : "",
                    "fileSizeInBytes" : ""
                    }],
               },
          }],
    }],
}

wgs_submission_grz = {
    "$schema": "https://raw.githubusercontent.com/BfArM-MVH/MVGenomseq/refs/tags/v1.2.1/GRZ/grz-schema.json",
    "submission" : {
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
        "labName" : ""
         },
    "donors" : [{
        "donorPseudonym" : "",
        "gender" : "",
        "relation" : "",
        "mvConsent" : {
            #"presentationDate" : "",
            "version" : "",
            "scope" : [{
                "type" : "",
                "date" : "",
                "domain" : "mvSequencing"
                },
                {
                "type" : "",
                "date" : "",
                "domain" : "reIdentification"
                },
                {
                "type" : "",
                "date" : "",
                "domain" : "caseIdentification"
                }
            ],
        },
        "researchConsents" : [
        ],
        #"researchConsents" : [{
        #    "schemaVersion" : "",
        #    "presentationDate" : "",
        #    "scope" : ""
        #    }],
        "labData" : [{
            "labDataName" : "",
            "tissueOntology" : {
                "name" : "",
                "version" : ""
                },
            "tissueTypeId" : "",
            "tissueTypeName" : "",
            "sampleDate" : "",
            "sampleConservation" : "",
            "sequenceType" : "",
            "sequenceSubtype" : "",
            "fragmentationMethod" : "",
            "libraryType" : "",
            "libraryPrepKit" : "",
            "libraryPrepKitManufacturer" : "",
            "sequencerModel" : "",
            "sequencerManufacturer" : "",
            "kitName" : "",
            "kitManufacturer" : "",
            "enrichmentKitManufacturer" : "",
            "enrichmentKitDescription" : "",
            "barcode" : "",
            "sequencingLayout" : "",
            "sequenceData" : {
                "bioinformaticsPipelineName" : "",
                "bioinformaticsPipelineVersion" : "",
                "referenceGenome" : "",
                "percentBasesAboveQualityThreshold" : {
                    "minimumQuality": "",
                    "percent": ""
                 },
                "meanDepthOfCoverage" : "",
                "minCoverage" : "",
                "targetedRegionsAboveMinCoverage" : "",
                "nonCodingVariants" : "",
                "callerUsed" : [{
                    "name" : "",
                   "version" : ""
                    }],
                "files" : [{
                    "filePath" : "",
                    "fileType" : "",
                    "checksumType" :"",
                    "fileChecksum" : "",
                    "fileSizeInBytes" : "",
                    "readOrder" : "",
                    "readLength" : "",
                    #"flowcellId" : "",
                    #"laneId" : ""
                    },
                    {
                    "filePath" : "",
                    "fileType" : "",
                    "checksumType" :"",
                    "fileChecksum" : "",
                    "fileSizeInBytes" : "",
                    "readOrder" : "",
                    "readLength" : "",
                    #"flowcellId" : "",
                    #"laneId" : ""
                    }],
               },
          },
          {
            "labDataName" : "",
            "tissueOntology" : {
                "name" : "",
                "version" : ""
                },
            "tissueTypeId" : "",
            "tissueTypeName" : "",
            "sampleDate" : "",
            "sampleConservation" : "",
            "sequenceType" : "",
            "sequenceSubtype" : "",
            "fragmentationMethod" : "",
            "libraryType" : "",
            "libraryPrepKit" : "",
            "libraryPrepKitManufacturer" : "",
            "sequencerModel" : "",
            "sequencerManufacturer" : "",
            "kitName" : "",
            "kitManufacturer" : "",
            "enrichmentKitManufacturer" : "",
            "enrichmentKitDescription" : "",
            "barcode" : "",
            "sequencingLayout" : "",
            "tumorCellCount" : [{
                "count" : "",
                "method" : ""
                }],
            "sequenceData" : {
                "bioinformaticsPipelineName" : "",
                "bioinformaticsPipelineVersion" : "",
                "referenceGenome" : "",
                "percentBasesAboveQualityThreshold" : {
                     "minimumQuality": "",
                     "percent": ""
                },
                "meanDepthOfCoverage" : "",
                "minCoverage" : "",
                "targetedRegionsAboveMinCoverage" : "",
                "nonCodingVariants" : "",
                "callerUsed" : [{
                    "name" : "",
                    "version" : ""
                    }],
                "files" : [{
                    "filePath" : "",
                    "fileType" : "",
                    "checksumType" :"",
                    "fileChecksum" : "",
                    "fileSizeInBytes" : "",
                    "readOrder" : "",
                    "readLength" : "",
                    #"flowcellId" : "",
                    #"laneId" : ""
                    },
                    {
                    "filePath" : "",
                    "fileType" : "",
                    "checksumType" :"",
                    "fileChecksum" : "",
                    "fileSizeInBytes" : "",
                    "readOrder" : "",
                    "readLength" : "",
                    #"flowcellId" : "",
                    #"laneId" : ""
                    },
                    {
                    "filePath" : "",
                    "fileType" : "",
                    "checksumType" :"",
                    "fileChecksum" : "",
                    "fileSizeInBytes" : ""
                    }],
               },
          }],
    }],
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
            #print(key)

            for fastq_file in lst_R1_R2_dict[index_read_fq][key].values():
                #print(fastq_file)

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

# login
def get_token(url):
    url = url
    headers = {
         "Content-Type":"application/json",
         "Authorization": gv.authorization
    }
    parameters = {}
    response = requests.post(url, headers=headers, json=parameters, verify=gv.certification)

    status = response.status_code
    response_json = response.json()
    token = response_json["response"]["token"]

    return token, status

# extract data
def get_all_records(token):
    url = gv.url_get_all_records
    headers = {"Authorization": "Bearer " + token}

    response = requests.get(url, headers=headers, verify=gv.certification)

    #print(response.json())

def get_record_Leistungserfassung(token,prefix,number,year,analysis):
    url = gv.url_get_record_leistungserfassung
    headers = {"Authorization": "Bearer "+ token}
    parameters = { "query":[
                     {  "Kürzel": prefix,
                        "Nummer": number,
                        "Jahr":year,
                        "Untersuchung":analysis}
                   ]
                 }


    response = requests.post(url, headers=headers, json=parameters, verify=gv.certification)

    #print(response.json())
    return response.json()

def get_record_nexus_mssql_db(token,OrderNumber):
    url = gv.url_get_record_nexus_mssql_db
    headers = {"Authorization": "Bearer "+ token}
    parameters = { "query":[
                     {  "OrderNumber": OrderNumber }
                   ]
                 }


    response = requests.post(url, headers=headers, json=parameters, verify=gv.certification)

    #print(response.json())

def get_orbis_id(token,lastname,firstname,renamed_dob):
    url = gv.url_get_orbis_id
    headers = {"Authorization": "Bearer "+ token}
    parameters = { "query":[
                     {  "Lastname": lastname,
                        "Firstname": firstname,
                         "dateofbirth": renamed_dob }
                   ]
                 }


    response = requests.post(url, headers=headers, json=parameters, verify=gv.certification)

    #print(response.json())
    return response.json()

def adjust_id_for_fmrest_request(sample_id):
    one_letter_prefix = ["J","M","P","G"] # E is excluded since it is not used for sample labeling - Exception!
    two_letter_prefix = ["JS","MS"]
    if sample_id[0] not in one_letter_prefix or sample_id[0] =="E":
        #prefix_ad_sample_id = "E-"+sample_id
        sample_id_splitted = sample_id.split("-")
        # since sequencing started not before 1999, the year 20?? is assumed
        sample_id_fmrest = ("E",sample_id_splitted[0],"20"+sample_id_splitted[1])

    elif sample_id[0] in one_letter_prefix and sample_id[0:2] not in two_letter_prefix:
        sample_id_splitted = sample_id[1:].split("-")
        sample_id_fmrest = (sample_id[0],sample_id_splitted[0],"20"+sample_id_splitted[1])

    elif sample_id[0:2] in two_letter_prefix:
        sample_id_splitted = sample_id[2:].split("-")
        sample_id_fmrest = (sample_id[0:2],sample_id_splitted[0],"20"+sample_id_splitted[1])

    return sample_id_fmrest

# logout

def parse_icdo3_xml(xml_file,code):
    tree = ET.parse(xml_file)
    root = tree.getroot()

    # check code
    if code.startswith("C") and len(code) <= 5:

        for child in root.findall("Class"):

            if code.count(".") == 1 and child.get("code") == code:
                if child[1].attrib["kind"] == "preferred":
                    for label in child[1]:
                        icdo3_text = label.text

            elif code.count(".") == 0 and child.get("code") == "C22":#code:
                if child[3].attrib["kind"] == "preferred":
                    for label in child[3]:
                        icdo3_text = label.text

    else:
        raise ValueError("ICD-O-3 code not valid!")

    return icdo3_text
