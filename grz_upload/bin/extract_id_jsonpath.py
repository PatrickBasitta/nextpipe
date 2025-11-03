#!/usr/bin/env python
import requests
import json
import functions_grz.grz_functions as grz
import functions_grz.grz_global_variables as gv
import os
import pandas as pd
import argparse
import re
from pathlib import Path
import glob
#import functions_grz.etch_research_consent as fsc

# using argparse for positinal arguments
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--target_dir_jsons", type=str)
parser.add_argument("-x", "--id_json_paths", type=str)
#parser.add_argument("-i", "--id_patient", type=str)
args = parser.parse_args()

# ETL path: json done
#parent_directory = os.path.abspath('.')
parent_directory = os.path.abspath(args.target_dir_jsons)
#target_dirs = [f for f in os.listdir(parent_directory) if re.match(r'[0-9]', f)]
files_DONE = glob.glob(parent_directory + "/**/" + "*_fertig.json", recursive=True)

print(files_DONE)

# make id_path_df for json generator
size_of_df = len(files_DONE)
id_path_df = pd.DataFrame("", index=list(range(0,size_of_df)),
                          columns = ["sample_id", 
                                     "json_path",
                                     "pid"])


for i, json_file in enumerate(files_DONE):
    print(i, json_file)
    
    # get case id and set xlsx_path
    basename = os.path.basename(json_file)
    sample_id = basename.split("_")[0] # get check !!!
    # check
    if basename.count(sample_id) == 1:
        id_path_df.loc[i,"sample_id"] = str(sample_id)
    else:
        raise ValueError(f"Please check the following patient '{sample_id}'")

    # get PID
    jfile=json_file
    with open(jfile, 'r') as file:
        data = json.load(file)
    print(data["patient_id"])

    id_path_df.loc[i,"json_path"] = Path(json_file)
    id_path_df.loc[i,"pid"] = data["patient_id"]
    print(id_path_df["sample_id"])
    print(id_path_df["json_path"])

# make samplesheet
id_json_df = id_path_df[["sample_id", "json_path", "pid"]]
id_json_df.to_csv(args.id_json_paths, sep=",",index=False)

