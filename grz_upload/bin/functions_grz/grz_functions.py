#!/usr/bin/env python

import requests
import functions_grz.grz_global_variables as gv
import xml.etree.ElementTree as ET
import subprocess
import sys

def get_broad_consent(url):
    url = url
    response = requests.get(url, auth=gv.authorization_bc, verify=gv.certification)

    return response.json()

def get_uuid(url,pid):
    url_get_uuid = url+pid
    response = requests.get(url_get_uuid, auth=gv.authorization_uuid, verify=gv.certification)

    return response.json()

def get_vn(vn_url_grz_knvr_uuid):

    response = requests.get(vn_url_grz_knvr_uuid, auth=gv.authorization_uuid, verify=gv.certification)

    return response.json()

def get_mv_consent(url_mv,pid):
    url_get_mv = url_mv+pid
    response = requests.get(url_get_mv, auth=gv.authorization_uuid, verify=gv.certification)

    return response.json()

def get_patient_pseudonym(url_patient_pseudonym,pid):
    suffix = gv.url_suffix_patient_pseudonym
    url_get_pp = url_patient_pseudonym+pid+suffix
    response = requests.get(url_get_pp, auth=gv.authorization_uuid, verify=gv.certification)
    text = response.text

    return text

def get_research_consent_pseudonym(url_research_consent_pseudonym,pid):
    suffix_bc = gv.url_suffix_research_consent_pseudonym
    url_get_bp = url_research_consent_pseudonym+pid+suffix_bc

    bc_pseudo_js = subprocess.run(
    [sys.executable, "fetch_research_consent.py",
                     "--patient-id", pid,
                     "--username", gv.username,
                     "--password", gv.password,
                     "--base-url", url_get_bp,
                     "--output", pid+"_bc.json",
                     "--ca-cert", gv.certification],
    capture_output = True,
    text = True,
    check = True
    )

def parse_icdo3_xml(xml_file,code):
    tree = ET.parse(xml_file)
    root = tree.getroot()

    code_lst = []
    for child in root.findall("Class"):
         code_lst.append(child.get("code"))

    # check code
    if code.startswith("C") and len(code) == 5:

        for child in root.findall("Class"):

            if code.count(".") == 1 and child.get("code") == code:
                print(child[1].attrib)
                if child[1].attrib["kind"] == "preferred":
                    for label in child[1]:
                        icdo3_text = label.text

            elif code not in code_lst:
                #print(child[3].attrib)
                #print(child[1].attrib)
                #if child[3].attrib["kind"] == "preferred":
                   #for label in child[3]:
                icdo3_text = "This code is not a ICD-O-3 code!!!"

    else:
        icdo3_text = "This code is not a ICD-O-3 code!!!"
        #raise ValueError("ICD-O-3 code not valid!")


    return icdo3_text
