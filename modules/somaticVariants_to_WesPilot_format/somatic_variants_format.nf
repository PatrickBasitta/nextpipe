nextflow.enable.dsl=2
/*
 * INPUT PARAMETERS
*/
params.input_file = ""
params.reference = ""
params.outdir = ""
params.case = ""
params.cpus = 12
params.memory = 60

/*
 * CHANNELS
*/
input_ch = Channel.fromPath(params.input_file)
chainhg38Tohg19_ch = Channel.fromPath(params.reference)

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
    processed_data = processed_data.rename(columns={"Chromosome": "Chr",
                                                "Position": "Start",
                                                "End Position": "End"})
    # write bed file
    bed_file = processed_data[["Chr", "Start", "End"]]
    bed_file.to_csv("output.bed", sep="\t", index=False)
    """
}

process FASTREMAP {
    
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
    case = "${params.case}"
    processed_data = pandas.read_excel(filename)
    bed_file = pandas.read_csv(bed, delimiter="\t", \
                               header=None, names=["Chr", "Start", "End"],index_col=False)
   
    # get hg19 coordinates and assing to processed data  
    if processed_data["Chromosome"].values in bed_file["Chr"].values:
        processed_data["Chr"] = bed_file["Chr"]
        processed_data["Start"] = bed_file["Start"]
        processed_data["End"] = bed_file["End"]

    # get WES Pilot format
    processed_data_final = processed_data[["Chr", "Start", "End",
                                     "Reference", "Allele", "HGVS_PROTEIN", \
                                     "GENE_SYMBOL", "TRANSCRIPT_ID", \
                                     "HGVS_TRANSCRIPT", "ING_AF", "Coverage_x"]]

    processed_data_final.insert(loc=0, column="Standort", value="UK_Bonn")
    processed_data_final.insert(loc=1, column="Sample", value=case)
    processed_data_final.insert(loc=11, column="CN", value="-")

    # write output
    processed_data_final.to_csv(case+".csv")
    """
}

workflow {
    MAKE_BED(input_ch)
    FASTREMAP(chainhg38Tohg19_ch,MAKE_BED.out.bed)
    WES_PILOT_FORMAT(input_ch,FASTREMAP.out.remapped_to_hg19)
    WES_PILOT_FORMAT.out.final_csv.view()
}


