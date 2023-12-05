nextflow.enable.dsl=2
/*
 * INPUT PARAMETERS
*/
params.input_file = ""
params.reference_chain = ""
params.segm_file= ""
params.outdir = ""
params.case_no = ""
params.cpus = 12
params.memory = 60

/*
 * CHANNELS
*/
input_ch = Channel.fromPath(params.input_file)
chainhg38Tohg19_ch = Channel.fromPath(params.reference_chain)

process MAKE_BED {
 
    input:
        path(input)
 
    output:
        path "output.bed", emit: bed
 
    """
    #!/usr/bin/env python
    import pandas

    filename = "${input}"

    processed_data = pandas.read_excel(filename)

    # rename Chromosome, Position, End Position to Chr, Start and End
    processed_data = processed_data.rename(columns={"Chromosome": "Chr", \
                                                "Position": "Start", \
                                                "End Position": "End"})
    # write bed file
    bed_file = processed_data[["Chr", "Start", "End"]]
    bed_file.to_csv("output.bed", sep="\t", index=False)
    """
}

process FASTREMAP {

    conda "bioconda::fastremap-bio" //Did not find corresponding container at quay.biocontainer
        
    cpus "${params.cpus}"
    memory "${params.memory}"
    debug true

    input:
       path(chainhg38Tohg19)
       path(bedfile)
 
    output:
       path("remapped_to_hg19.bed"), emit: remapped_to_hg19

    """
    FastRemap -f bed -c ${chainhg38Tohg19} -i ${bedfile} -u unmapped.bed -o remapped_to_hg19
    """
} 

process WES_PILOT_FORMAT {
    
    publishDir "${params.outdir}", mode: "copy"

    input:
        path(input)
        path(remapped_tohg19_bedfile)
    
    output:
        path "*.csv", emit: final_csv
   
     """
    #!/usr/bin/env python
    import pandas

    filename = "${input}"
    bed = "${remapped_tohg19_bedfile}"
    case = "${params.case_no}"
    segment_file = "${params.segm_file}"
    processed_data = pandas.read_excel(filename)
    bed_file = pandas.read_csv(bed, delimiter="\t", \
                               header=None, names=["Chr", "Start", "End"],index_col=False)
   
    # get hg19 coordinates and assing to processed data  
    if processed_data["Chromosome"].values in bed_file["Chr"].values:
        processed_data["Chr"] = bed_file["Chr"]
        processed_data["Start"] = bed_file["Start"]
        processed_data["End"] = bed_file["End"]

    # get WES Pilot format
    processed_data_final = processed_data[["Chromosome", "Position", "End Position", \
                                     "Reference", "Allele", "GENE_SYMBOL", "HGVS_PROTEIN", \
                                     "TRANSCRIPT_ID", \
                                     "HGVS_TRANSCRIPT", "Frequency_x", "Coverage_x"]]
    
     # round AF to max 2 decimals (since AF number are str in excel, one has to set these values to int before run this script!
    processed_data_final.loc[:,"Frequency_x"]  = processed_data_final\
                                   ["Frequency_x"].apply(lambda x: round(x,2))

    # rename according to WES_Pilot_format
    processed_data_final = processed_data_final.rename(columns={"Chromosome": "Chr", \
                                                                "Position": "Start", \
                                                                "End Position": "End", \
                                                                "Reference": "Ref", \
                                                                "Allele": "Alt", \
                                                                "GENE_SYMBOL": "Gen", \
                                                                "HGVS_PROTEIN": "Aminosaeure", \
                                                                "TRANSCRIPT_ID": "NM-Nummer", \
                                                                "HGVS_TRANSCRIPT": "cDNA", \
                                                                "Frequency_x": "VAF", \
                                                                "Coverage_x": "Coverage"})    

    processed_data_final.insert(loc=0, column="Standort", value="Bonn")
    processed_data_final.insert(loc=1, column="Sample", value=case)
    processed_data_final.insert(loc=11, column="CN", value="-")

    # get CN number from segment_file
    segm_data = pandas.read_csv(segment_file, delimiter= "\t", engine="python")

    # rename
    segm_data = segm_data.rename(columns={"chromosome": "Chr", \
                                          "start.pos": "Start", \
                                          "end.pos": "End"})

    # get chr from variants
    for i in range(len(processed_data_final)):
    
        # get corresponding chr
        segm_data_chr = segm_data[segm_data["Chr"]==str(processed_data_final["Chr"][i])]

        # get df with CN value
        tmp_df = processed_data_final.iloc[i,:].to_frame().T
    
        df_check_range = tmp_df.merge\
                        (segm_data_chr, how="cross").query\
                        ("(Start_x >= Start_y) & (End_x <= End_y)")
    
   
        if df_check_range.shape[0] == 1:
            cn_value = df_check_range["CNt"][df_check_range.index[0]]
            processed_data_final.loc[i,"CN"] = cn_value

    # write output
    processed_data_final.to_csv(case+".csv", index=False)
    """
}

workflow {
    MAKE_BED(input_ch)
    FASTREMAP(chainhg38Tohg19_ch,MAKE_BED.out.bed)
    WES_PILOT_FORMAT(input_ch,FASTREMAP.out.remapped_to_hg19)
    WES_PILOT_FORMAT.out.final_csv.view()
}


