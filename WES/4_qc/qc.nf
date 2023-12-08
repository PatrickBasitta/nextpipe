nextflow.enable.dsl=2
/*
 * prelimnary script to polish CLC bam files (exlude cigar string error) and 
 * perform picard CollectHsMetrics for WES Pilot
 */

/*
 * Input files
 */
params.bamfile = "/PAT-Sequenzer/ImportExport/BAM_WES/test/*_clc_*.bam"
//params.qc_csv = ""
params.outdir = ""
params.CLC38_fa = ""
params.TB_list = ""
params.cpus = 12
params.memory = 60

/*
 * Channels
 */

//qc_csv_ch = Channel.fromPath(params.qc_csv) | splitCsv(header: true) | map { row -> tuple(row.sampleid,row.bamfile) }
samples_ch = Channel.fromPath(params.bamfile) | map { tuple(it) }
CLC38_fa_ch = Channel.fromPath(params.CLC38_fa)
TB_list_ch = Channel.fromPath(params.TB_list)

process POLISH_BAM {
    
    tag "${bamfile.getSimpleName()}"
    cpus "${params.cpus}"
    memory "${params.memory}"     
    debug true
    //publishDir "${params.outdir}/QC/${bamfile.getSimpleName()}", mode: "copy"      
    input:
        path(bamfile)
    output:
        path("${bamfile.getSimpleName()}_polished.bam"), emit: bamout
    script:
    """
    samtools view -h ${bamfile} | grep -v ".S.P" | grep -v ".S..P" | grep -v ".S...P" | samtools view -S -b -o ${bamfile.getSimpleName()}_polished.bam - 
    """
} 

process PICARD_CHM {
    
    tag "${bamout.getSimpleName()}"
    cpus "${params.cpus}"
    memory "${params.memory}"
    debug true
    publishDir "${params.outdir}/QC/${bamout.getSimpleName()}", mode: "copy"
    input:
        path(bamout)
        path(CLC38_fa)
        path(TB_list)
    output:
        path("${bamout.getSimpleName()}.csv")
    script:
    """
    mkdir -p QC/${bamout.getSimpleName()}
    picard CollectHsMetrics -I ${bamout} -O ${bamout.getSimpleName()}.csv -R ${CLC38_fa} -TARGET_INTERVALS ${TB_list} -BAIT_INTERVALS ${TB_list} 
    """
}

workflow QC_WES_PILOT {
    take:
        samples_ch
        CLC38_fa_ch
        TB_list_ch  
    main:
        POLISH_BAM(samples_ch)
        PICARD_CHM(POLISH_BAM.out.bamout,CLC38_fa_ch,TB_list_ch)
    emit:
        PICARD_CHM.out
}
    
workflow {
    //samples_ch.view()
    QC_WES_PILOT(samples_ch,CLC38_fa_ch,TB_list_ch)
    //QC_WES_PILOT.out.view()
}
