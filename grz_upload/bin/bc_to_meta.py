#!/usr/bin/env python
import requests
import json
import functions_grz.grz_functions as grz
import functions_grz.grz_global_variables as gv
import argparse
import xml.etree.ElementTree as ET

# using argparse for positinal arguments
parser = argparse.ArgumentParser()
parser.add_argument("-b", "--bc_file", type=str)
parser.add_argument("-m", "--meta_json", type=str)
parser.add_argument("-i", "--sample_id", type=str)
parser.add_argument("-s", "--genomic_study_subtype", type=str)
parser.add_argument("-l", "--library_type", type=str)
parser.add_argument("-p", "--patient_id", type=str)
parser.add_argument("-x", "--icdo3_xml", type=str)
args = parser.parse_args()

# load data
with open(args.bc_file, 'r') as bcf:
    bc_data = json.load(bcf)

with open(args.meta_json, 'r') as mj:
    mj_data = json.load(mj)

# pid
pid = str(args.patient_id)

# add mv consent
mv_consent = grz.get_mv_consent(gv.url_mv,pid)
mv_consent_data_wo_presentation_date = {j: i for j, i in mv_consent.items() if j != "presentationDate"}
mj_data["submission"]["donors"][0]["mvConsent"] = mv_consent_data_wo_presentation_date

# add bc consent
mj_data["submission"]["donors"][0]["researchConsents"] = [bc_data]

# removes metaData UKB and KDK part
#mj_data["analysis"]
#mj_data["pathoProId"]
#mj_data["nexusId"]
#mj_data["molpathId"]
#mj_data["PID"]
#mj_data["firstname"]
#mj_data["lastame"]
#mj_data["dateofbrith"]
#mj_data["network"]
#mj_data["diseaseType"]
#mj_data["entity_excel"]
#mj_data["entity_fm"]
#mj_data["molecularReportod"]

# removes labData[0] == normal part
#if args.genomic_study_subtype == "tumor-only" and args.library_type == "wes":
#    del mj_data["submission"]["donors"][0]["labData"][0]
#print(json.dumps(bc_data,indent=4))
#print(json.dumps(mj_data["submission"]["donors"][0],indent=4))
#print(json.dumps(mj_data,indent=4))

# get uuid
#pid = str(args.patient_id)
#print(pid)
uuid_url = gv.url_uuid
#print(uuid_url)
#url = str(uuid_url+pid)
#print(url)
uuid_block = grz.get_uuid(uuid_url,pid)
print(uuid_block)
#print(json.dumps(uuid_block,indent=4))
#print(uuid_block)
if pid == str(uuid_block[0]["patient_id"]):
    uuid = uuid_block[0]["process_id"]
else:
    raise ValueError("PIDs do not match!!!")

# add uuid
mj_data["mvtracker_uuid"] = uuid
mj_data["submission"]["submission"]["localCaseId"] = uuid
#mj_data["submission"]["donors"][0]["donorPseudonym"] = uuid

# get tanG
if pid == str(uuid_block[0]["patient_id"]):
    vn = uuid_block[0]["kvnr"]
else:
    raise ValueError("PIDs do not match!!!")

#grz_data_center
vn_url = gv.url_vn
vn_url_grz_knvr_uuid = vn_url+vn+"&an="+uuid
vn_block = grz.get_vn(vn_url_grz_knvr_uuid)

if vn_block["success"]:
    tang = vn_block["vn"]
else:
    raise valueError("response error")

mj_data["submission"]["submission"]["tanG"] = tang

# add icd-o-3
if args.genomic_study_subtype == "tumor+germline" or args.genomic_study_subtype == "tumor-only" and args.library_type == "wxs":
    tumor_index = 1
    xml_file = args.icdo3_xml
    #print(xml_file)
    code = mj_data["submission"]["donors"][0]["labData"][tumor_index]["tissueTypeId"]
    #print(code)
    icdo3_txt = grz.parse_icdo3_xml(xml_file,code)
    mj_data["submission"]["donors"][0]["labData"][tumor_index]["tissueTypeName"] = icdo3_txt

if args.genomic_study_subtype == "tumor-only" and args.library_type == "panel":
    xml_file = args.icdo3_xml
    #print(xml_file)
    code = mj_data["submission"]["donors"][0]["labData"][0]["tissueTypeId"]
    #print(code)
    icdo3_txt = grz.parse_icdo3_xml(xml_file,code)
    mj_data["submission"]["donors"][0]["labData"][0]["tissueTypeName"] = icdo3_txt

#if args.genomic_study_subtype == "tumor-only" and args.library_type == "wxs":
#    tumor_index = 1
#    xml_file = args.icdo3_xml
    #print(xml_file)
#    code = mj_data["submission"]["donors"][0]["labData"][tumor_index]["tissueTypeId"]
    #print(code)
#    icdo3_txt = grz.parse_icdo3_xml(xml_file,code)
#    mj_data["submission"]["donors"][0]["labData"][tumor_index]["tissueTypeName"] = icdo3_txt

# patho meta - all
patho_meta = mj_data
with open(args.sample_id + "_patho_meta.json", "w") as file:
    json.dump(
    patho_meta,
    file,
    indent=4,
    sort_keys=False
)

# grz meta - submission
# removes labData[0] == normal part
# removes vcf file at index file[2]
if args.genomic_study_subtype == "tumor-only" and args.library_type == "wxs":
    del mj_data["submission"]["donors"][0]["labData"][0]
    # tumor index 0 in tumor-only mode - remove vcf
    del mj_data["submission"]["donors"][0]["labData"][0]["sequenceData"]["files"][2]

# in normal-tumor mode - remove vcf
if args.genomic_study_subtype == "tumor+germline" and args.library_type == "wxs":
    #normal_index = 0
    tumor_index =1
    #del mj_data["submission"]["donors"][0]["labData"][normal_index]["sequenceData"]["files"][2]
    del mj_data["submission"]["donors"][0]["labData"][tumor_index]["sequenceData"]["files"][2]

# in tumor-only and panel mode - remove vcf
if args.genomic_study_subtype == "tumor-only" and args.library_type == "panel":
    del mj_data["submission"]["donors"][0]["labData"][0]["sequenceData"]["files"][2]

# grz submission grz
grz_meta = mj_data["submission"]
with open("metadata.json", "w") as file:
    json.dump(
    grz_meta,
    file,
    indent=4,
    sort_keys=False
)

