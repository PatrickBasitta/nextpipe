nextflow.enable.dsl=2
/*
 * To run this nxf workflow the ref_fasta has to bgziped first using bzip -c .fasta > fasta.gz
 * and the command sequenza-utils gc_wiggle has to be applied to create .wig.gz file as input
 */

params.cpus = 12
params.memory = 60
params.outdir = ""
params.input_csv = ""
params.ref_fasta_gz = ""
params.gc_wig_file = ""

process BAM2SEQZ {
    debug true     
    cpus "${params.cpus}"
    memory "${params.memory}"

    publishDir "${params.outdir}/${id}/2_bam2seqz", mode: "copy"

    input:
        tuple val(id), file(normal_bam), file(tumor_bam) 
        path(ref_fasta_gz)
        path(gc_wig_file)

    output:
       tuple val(id) path("${id}"_out_seqz.gz)

    script:
    """
    mkdir -p ${id}/2_bam2seqz
    sequenza-utils bam2seqz --normal ${normal_bam} --tumor ${tumor_bam} --fasta ${ref_fasta_gz} -gc ${gc_wig_file} --output ${id}_out.seqz.gz
    """
}

process SEQZ_BINNING {

    cpus "${params.cpus}"
    memory "${params.memory}"

    publishDir "${params.outdir}/${id}/2_bam2seqz", mode: "copy"

    input:
        tuple val(id) path(bam2seqfile)
        
    output:
        tuple val(id) path("${id}"_bin50.seqz.gz)

    script:
    """
    sequenza-utils seqz_binning --seqz ${bam2seqfile} --window 50 -o ${id}_bin50.seqz.gz
    """
}

process R_SEQUENZA { 
/*
 * manually installed conda package r-sequenza and corresponding packages via mamba installer and BiocManager 
 * (here GenomeinfoDbData) locally
 */
    cpus "${params.cpus}"
    memory "${params.memory}"
    publishDir "${params.outdir}/${id}", mode: "copy"
    
    input:
        tuple val(id) path(binfile)
  
    output:
       path("*")

    script:
    """
    #!/usr/bin/env Rscript
    require("sequenza")
    Sys.setenv("VROOM_CONNECTION_SIZE" = 5000072 * 250)
    data.file <- ("${binfile}")
    seqzdata <- sequenza.extract(data.file)
    CP <- sequenza.fit(seqzdata)
    dir.create("3_final_results")
    sequenza.results(sequenza.extract = seqzdata, cp.table = CP, sample.id = "${id}", out.dir = "3_final_results")
    """
}     

/*
 * SUBWORKFLOW
 */

workflow sequenza_CLC {
      
    take:
        csv_ch
        ref_fasta_gz_ch
        gc_wig_file_ch
                
    main:
       BAM2SEQZ(csv_ch,ref_fasta_gz_ch,gc_wig_file_ch) \
       SEQZ_BINNING(BAM2SEQZ.out) \
       R_SEQUENZA(SEQZ_BINNING.out)
}

/*
 * MAIN WORKFLOW
 */

 csv_ch = Channel.fromPath(params.input_csv) | splitCsv(header: true) | map { row-> tuple(row.sampleid, file(row.normal_bam), file(row.tumor_bam)) }
 ref_fasta_gz_ch = Channel.fromPath(params.ref_fasta_gz)
 gc_wig_file_ch = Channel.fromPath(params.gc_wig_file)

workflow {   
    take:
        csv_ch
        ref_fasta_gz_ch
        gc_wig_file_ch

    main:
        sequenza_CLC(csv_ch,ref_fasta_gz_ch,gc_wig_file_ch)
}
