nextflow.enable.dsl=2
/*
 * To run this nxf workflow the ref_fasta has to bgziped first using bzip -c .fasta > fasta.gz
 * and the command sequenza-utils gc_wiggle has to be applied to create .wig.gz file as input
 */

/*
 * pipeline input parameters
 */
//params.case = ""
params.input_csv = ""
params.ref_fasta_gz = ""
params.gc_wig_file = ""
//params.normal_bam = ""
//params.tumor_bam = ""
params.outdir = ""
params.cpus = 12
params.memory = 60

/*
 * CHANNELS
 */

csv_ch_n = Channel.fromPath(params.input_csv) | splitCsv(header: true) | map { row -> tuple(row.sampleid,row.normalbam) }
csv_ch_t = Channel.fromPath(params.input_csv) | splitCsv(header: true) | map { row -> tuple(row.sampleid,row.tumorbam) }
//normal_bam_ch = Channel.fromPath(params.normal_bam)
//tumor_bam_ch = Channel.fromPath(params.tumor_bam)
ref_fasta_gz_ch = Channel.fromPath(params.ref_fasta_gz)
gc_wig_file_ch = Channel.fromPath(params.gc_wig_file)

process SAM_SORT_normal_BAM {
   
    tag "${sampleid}"    
    cpus "${params.cpus}"
    memory "${params.memory}"
    
    publishDir "${params.outdir}/${sampleid}/1_bai", mode: "copy"

    input:
        tuple val(sampleid), path(normalbam)

    output:
        tuple val(sampleid), path("*.bam")

    script:
    """
    mkdir -p ${sampleid}/1_bai
    samtools sort ${normalbam} -o ${sampleid}_sorted_normal.bam
    """
}

process normal_BAM_INDEX {
   
    cpus "${params.cpus}"
    memory "${params.memory}"

    publishDir "${params.outdir}/${sampleid}/1_bai", mode: "copy"

    input:
       tuple val(sampleid), path(sorted_n_bam)
       
    output:
        tuple val(sampleid), path("*.bai"), emit: normalbai

    script:
    """
    samtools index ${sorted_n_bam} ${sorted_n_bam}.bai 
    """
}

process SAM_SORT_tumor_BAM {
     
    tag "${sampleid}"
    cpus "${params.cpus}"
    memory "${params.memory}"

    publishDir "${params.outdir}/${sampleid}/1_bai", mode: "copy"

    input:
        tuple val(sampleid), path(tumorbam)

    output:
        tuple val(sampleid), path("*.bam")

    script:
    """
    mkdir -p ${sampleid}/1_bai
    samtools sort ${tumorbam} -o ${sampleid}_sorted_tumor.bam
    """
}

process tumor_BAM_INDEX {

    cpus "${params.cpus}"
    memory "${params.memory}"

    publishDir "${params.outdir}/${sampleid}/1_bai", mode: "copy"

    input:
        tuple val(sampleid), path(sorted_t_bam)

    output:
        tuple val(sampleid), path("*.bai"), emit: tumorbai

    script:
    """
    samtools index ${sorted_t_bam} ${sorted_t_bam}.bai
    """
}

process BAM2SEQZ {

    cpus "${params.cpus}"
    memory "${params.memory}"

    publishDir "${params.outdir}/${sampleid}/2_bam2seqz", mode: "copy"

    input:
        tuple val(sampleid), path(sorted_n_bam)
        tuple val(sampleid), path(sorted_t_bam)
        path(ref_fasta_gz)
        path(gc_wig_file)

    output:
       tuple val(sampleid), path("*.gz")

    script:
    """
    mkdir -p ${sampleid}/2_bam2seqz
    sequenza-utils bam2seqz --normal ${sorted_n_bam} --tumor ${sorted_t_bam} --fasta ${ref_fasta_gz} -gc ${gc_wig_file} --output ${sampleid}_out.seqz.gz
    """
}

process SEQZ_BINNING {

    cpus "${params.cpus}"
    memory "${params.memory}"

    publishDir "${params.outdir}/${sampleid}/2_bam2seqz", mode: "copy"

    input:
        tuple val(sampleid), path(bam2seqfile)
        
    output:
        tuple val(sampleid), path("*seqz.gz")

    script:
    """
    sequenza-utils seqz_binning --seqz ${bam2seqfile} --window 50 -o ${sampleid}_bin50.seqz.gz
    """
}

process R_SEQUENZA { 
/*
 * manually installed conda package r-sequenza and corresponding packages via mamba installer and BiocManager 
 * (here GenomeinfoDbData) locally
 */
    cpus "${params.cpus}"
    memory "${params.memory}"
    publishDir "${params.outdir}/${sampleid}", mode: "copy"
    
    input:
        tuple val(sampleid), path(binfile)
  
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
    sequenza.results(sequenza.extract = seqzdata, cp.table = CP, sample.id = "${sampleid}", out.dir = "3_final_results")
    """
}     

workflow SEQUENZA {
    take:
        csv_ch_n
        csv_ch_t
        ref_fasta_gz_ch
        gc_wig_file_ch
    main:
    SAM_SORT_normal_BAM(csv_ch_n)
    SAM_SORT_tumor_BAM(csv_ch_t)
    normal_BAM_INDEX(SAM_SORT_normal_BAM.out)
    tumor_BAM_INDEX(SAM_SORT_tumor_BAM.out)
    //BAM2SEQZ(SAM_SORT_normal_BAM.out.join(SAM_SORT_tumor_BAM.out),ref_fasta_gz_ch,gc_wig_file_ch)
    BAM2SEQZ(SAM_SORT_normal_BAM.out, SAM_SORT_tumor_BAM.out,ref_fasta_gz_ch,gc_wig_file_ch)
    SEQZ_BINNING(BAM2SEQZ.out)
    R_SEQUENZA(SEQZ_BINNING.out)
}

workflow {
    SEQUENZA(csv_ch_n,csv_ch_t,ref_fasta_gz_ch,gc_wig_file_ch)
}
