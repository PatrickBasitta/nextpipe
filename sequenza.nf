nextflow.enable.dsl=2
/*
 * pipeline input parameters
 */
params.case = ""
params.normal_cram = ""
params.tumor_cram =""
params.ref_fasta = ""
params.ref_fasta_gz = ""
params.gc_wig_file = ""
params.outdir = ""
params.cpus = 12
params.memory = 60

/*
 * CHANNELS
 */

normal_cram_ch = Channel.fromPath(params.normal_cram)
tumor_cram_ch = Channel.fromPath(params.tumor_cram)
ref_fasta_ch = Channel.fromPath(params.ref_fasta)
ref_fasta_gz_ch = Channel.fromPath(params.ref_fasta_gz)
gc_wig_file_ch = Channel.fromPath(params.gc_wig_file)

process normal_CRAM_TO_BAM {
    
    cpus "${params.cpus}"
    memory "${params.memory}"
    
    publishDir "${params.outdir}/${params.case}/1_bam", mode: "copy"

    input:
        path(n_cram)
        path(ref_fasta)

    output:
        path("*.bam"), emit: normalbam 
        
    script:
    """
    mkdir -p ${params.case}/1_bam
    samtools view \
        -T ${ref_fasta} \
        -o ${params.case}_normal.bam \
        ${n_cram} 
    """
}

process tumor_CRAM_TO_BAM {

    cpus "${params.cpus}"
    memory "${params.memory}"

    publishDir "${params.outdir}/${params.case}/1_bam", mode: "copy"

    input:
        path(t_cram)
        path(ref_fasta)

    output:
        path("*.bam"), emit: tumorbam

    script:
    """
    samtools view \
        -T ${ref_fasta} \
        -o ${params.case}_tumor.bam \
        ${t_cram}
    """
}

process normal_BAM_INDEX {
   
    cpus "${params.cpus}"
    memory "${params.memory}"

    publishDir "${params.outdir}/${params.case}/1_bam", mode: "copy"

    input:
        path(n_bam)
       
    output:
        path("*.bai"), emit: normalbai

    script:
    """
    samtools index ${n_bam} ${n_bam}.bai 
    """
}

process tumor_BAM_INDEX {

    cpus "${params.cpus}"
    memory "${params.memory}"

    publishDir "${params.outdir}/${params.case}/1_bam", mode: "copy"

    input:
        path(t_bam)

    output:
        path("*.bai"), emit: tumorbai

    script:
    """
    samtools index ${t_bam} ${t_bam}.bai
    """
}

process BAM2SEQZ {

    cpus "${params.cpus}"
    memory "${params.memory}"

    publishDir "${params.outdir}/${params.case}/2_bam2seqz", mode: "copy"

    input:
        path(n_bam)
        path(t_bam)
        path(ref_fasta_gz)
        path(gc_wig_file)

    output:
        path("*.gz"), emit: bam2seqfile

    script:
    """
    mkdir -p ${params.case}/2_bam2seqz
    sequenza-utils bam2seqz --normal ${n_bam} --tumor ${t_bam} --fasta ${ref_fasta_gz} -gc ${gc_wig_file} --output ${params.case}_out.seqz.gz
    """
}

process SEQZ_BINNING {

    cpus "${params.cpus}"
    memory "${params.memory}"

    publishDir "${params.outdir}/${params.case}/2_bam2seqz", mode: "copy"

    input:
        path(bam2seqfile)
        
    output:
        path("*seqz.gz"), emit: binfile

    script:
    """
    sequenza-utils seqz_binning --seqz ${bam2seqfile} --window 50 -o ${params.case}_bin50.seqz.gz
    """
}

process R_SEQUENZA {
    
    container "docker.io/sequenza/sequenza:latest"
    
    cpus "${params.cpus}"
    memory "${params.memory}"
    stageInMode "copy"
    publishDir "${params.outdir}/${params.case}", mode: "copy"
    
    input:
        path(binfile)
  
    output:
       path("*")

    script:
    """
    #!/usr/bin/env Rscript
    library("sequenza")
    data.file <- "${binfile}"
    seqzdata <- sequenza.extract(data.file)
    CP <- sequenza.fit(seqzdata)
    sequenza.results(sequenza.extract = seqzdata, cp.table = CP, sample.id = "${params.case}", out.dir="/home/sequenza/")
    """
}     

workflow {
    normal_CRAM_TO_BAM(normal_cram_ch, ref_fasta_ch)  
    tumor_CRAM_TO_BAM(tumor_cram_ch, ref_fasta_ch)
    normal_BAM_INDEX(normal_CRAM_TO_BAM.out.normalbam)
    tumor_BAM_INDEX(tumor_CRAM_TO_BAM.out.tumorbam)
    BAM2SEQZ(normal_CRAM_TO_BAM.out.normalbam,tumor_CRAM_TO_BAM.out.tumorbam,ref_fasta_gz_ch,gc_wig_file_ch)
    SEQZ_BINNING(BAM2SEQZ.out.bam2seqfile)
    R_SEQUENZA(SEQZ_BINNING.out.binfile)
}

