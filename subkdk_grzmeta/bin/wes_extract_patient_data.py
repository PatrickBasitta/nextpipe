#!/usr/bin/env python

import argparse
import requests
import json as js
import functions_etl.MVdataset_generator_utils as etl
import functions_etl.global_variables as gv

# set arguments
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--patient_id", type=str)
parser.add_argument("-o", "--patient_meta_data", type=str)
args = parser.parse_args()

# Get token for FileMaker API Operations
url = gv.url_token
token,status = etl.get_token(url)

# Get record firstname, lastname, dob, orbis_id
# possible prefixes:
#    without prefix = E
#    others are: J,JS,M,MS,P,G
sample_id = args.patient_id
sample_id_fmrest = etl.adjust_id_for_fmrest_request(sample_id)
prefix = list(sample_id_fmrest)[0]
number = list(sample_id_fmrest)[1]
year = list(sample_id_fmrest)[2]
#print(prefix,number,year)
#print(sample_id_fmrest)
if status == 200:
    #print(token,status)
    #get_all_records(token)
    get_record_LE = etl.get_record_Leistungserfassung(token,prefix,number,year,"WES")
    #print(get_record_LE)
    if get_record_LE["response"]["dataInfo"]["foundCount"] == get_record_LE["response"]["dataInfo"]["returnedCount"]:
        #print(get_record_LE["response"]["data"][0]["fieldData"])
        firstname = get_record_LE["response"]["data"][0]["fieldData"]["Vorname"]
        lastname = get_record_LE["response"]["data"][0]["fieldData"]["Name"]
        dob = get_record_LE["response"]["data"][0]["fieldData"]["GBD"]
        analysis = get_record_LE["response"]["data"][0]["fieldData"]["Untersuchung"]
        entity = get_record_LE["response"]["data"][0]["fieldData"]["Diagnose"]
        print(firstname, lastname, dob)
        split_dob = dob.split("/")
        #print(split_dob)
        renamed_dob = split_dob[2]+split_dob[0]+split_dob[1]
        #print(renamed_dob)
        nexus_mssqldb_data = etl.get_orbis_id(token,lastname,firstname,renamed_dob) # plus id for better identification
        #print(nexus_mssqldb_data)
        orbis_id_lst = []
        gender_lst = []
        if nexus_mssqldb_data["response"]["dataInfo"]["foundCount"] == nexus_mssqldb_data["response"]["dataInfo"]["returnedCount"]:
            record_count = nexus_mssqldb_data["response"]["dataInfo"]["foundCount"]
            for i in range(record_count):
                #print(i)
                #print(nexus_mssqldb_data["response"]["data"][i]["fieldData"])
                patient_orbis_id = nexus_mssqldb_data["response"]["data"][i]["fieldData"]["Patienten ID"]
                orbis_id_lst.append(patient_orbis_id)
                # get gender
                gender = nexus_mssqldb_data["response"]["data"][i]["fieldData"]["Gender"]
                gender_lst.append(gender)

        if len(set(orbis_id_lst)) == 1:
            orbis_id = set(orbis_id_lst)
        else:
            #raise ValueError("Multiple orbis_ids available!")
            orbis_id = "?"
        #print(list(orbis_id)[0])

        if len(set(gender_lst)) == 1:
            if gender == 1:
                sex = "male"
            elif gender == 2:
                sex = "female"
            else:
                sex = "unknown"
        else:
            #raise ValueError("Gender is not clear, different entries exist!!!")
            sex = "Please check, because different entries exists!!!"
        # be aware: the value "other" is here excluded
        # print(sex)
    else:
        analysis = "WES"
        firstname = ""
        lastname = ""
        dob = ""
        orbis_id = "?"
        sex = ""
        entity = ""

    patient_data = {
                          "analysis": analysis,
                          "firstname": firstname,
                          "lastname": lastname,
                          "dateofbirth": dob,
                          "pathoProId": "/".join(sample_id_fmrest),
                          "nexusId": "",
                          "molpathId": "",
                          "orbisId": list(orbis_id)[0].replace(" ", ""),
                          "gender": sex,
                          "entity": entity
                          }

    with open(args.patient_meta_data, "w") as file:
        js.dump(
        patient_data,
        file,
        indent=4,
        sort_keys=False
        )
