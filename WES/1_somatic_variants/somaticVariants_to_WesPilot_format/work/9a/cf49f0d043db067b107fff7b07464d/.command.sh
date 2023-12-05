#!/usr/bin/env python
import pandas

filename = "intermediate_3_ready.xlsx"
bed = "remapped_to_hg19.bed"
case = "3"
segment_file = "/PAT-Sequenzer/WES-Pilot/WES_Pilot_input_all/CNV_all/3_segments_hg19.txt"
processed_data = pandas.read_excel(filename)
bed_file = pandas.read_csv(bed, delimiter="	",                                header=None, names=["Chr", "Start", "End"],index_col=False)

# get hg19 coordinates and assing to processed data  
if processed_data["Chromosome"].values in bed_file["Chr"].values:
    processed_data["Chr"] = bed_file["Chr"]
    processed_data["Start"] = bed_file["Start"]
    processed_data["End"] = bed_file["End"]

# get WES Pilot format
processed_data_final = processed_data[["Chromosome", "Position", "End Position",                                      "Reference", "Allele", "GENE_SYMBOL", "HGVS_PROTEIN",                                      "TRANSCRIPT_ID",                                      "HGVS_TRANSCRIPT", "Frequency_x", "Coverage_x"]]

 # round AF to max 2 decimals (since AF number are str in excel, one has to set these values to int before run this script!
processed_data_final.loc[:,"Frequency_x"]  = processed_data_final                                   ["Frequency_x"].apply(lambda x: round(x,2))

# rename according to WES_Pilot_format
processed_data_final = processed_data_final.rename(columns={"Chromosome": "Chr",                                                                 "Position": "Start",                                                                 "End Position": "End",                                                                 "Reference": "Ref",                                                                 "Allele": "Alt",                                                                 "GENE_SYMBOL": "Gen",                                                                 "HGVS_PROTEIN": "Aminosaeure",                                                                 "TRANSCRIPT_ID": "NM-Nummer",                                                                 "HGVS_TRANSCRIPT": "cDNA",                                                                 "Frequency_x": "VAF",                                                                 "Coverage_x": "Coverage"})    

processed_data_final.insert(loc=0, column="Standort", value="Bonn")
processed_data_final.insert(loc=1, column="Sample", value=case)
processed_data_final.insert(loc=11, column="CN", value="-")

# get CN number from segment_file
segm_data = pandas.read_csv(segment_file, delimiter= "	", engine="python")

# rename
segm_data = segm_data.rename(columns={"chromosome": "Chr",                                           "start.pos": "Start",                                           "end.pos": "End"})

# get chr from variants
for i in range(len(processed_data_final)):

    # get corresponding chr
    segm_data_chr = segm_data[segm_data["Chr"]==str(processed_data_final["Chr"][i])]

    # get df with CN value
    tmp_df = processed_data_final.iloc[i,:].to_frame().T

    df_check_range = tmp_df.merge                        (segm_data_chr, how="cross").query                        ("(Start_x >= Start_y) & (Ende_x <= Ende_y)")


    if df_check_range.shape[0] == 1:
        cn_value = df_check_range["CNt"][df_check_range.index[0]]
        processed_data_final.loc[i,"CN"] = cn_value

# write output
processed_data_final.to_csv(case+".csv", index=False)
