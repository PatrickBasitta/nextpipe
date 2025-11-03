#!/usr/bin/env python

import requests
import functions_grz.grz_global_variables as gv
import xml.etree.ElementTree as ET

def get_broad_consent(url):
    url = url
    response = requests.get(url, auth=gv.authorization_bc, verify=gv.certification)

    return response.json()

def get_uuid(url,pid):
    url_get_uuid = url+pid
    response = requests.get(url_get_uuid, auth=gv.authorization_uuid, verify=gv.certification)

    return response.json()

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

            elif code.count(".") == 0 and child.get("code") == code:
                if child[3].attrib["kind"] == "preferred":
                    for label in child[3]:
                        icdo3_text = label.text

            else:
                icdo3_text = "This code is not a ICD-O-3 code!!!"

    else:
        raise ValueError("ICD-O-3 code not valid!")

    return icdo3_text
