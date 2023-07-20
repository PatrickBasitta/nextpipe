nextflow.enable.dsl=2
/*
 * pipeline input parameters
 */
params.reads = "/data/data_AG_Glen/input_data/DE65NGSUKBD135763_3-LNCaP-Apa_{1,2}.fq.gz"
params.reference = "/data/data_AG_Glen/reference_data/"
params.outdir = "results"
params.cpus = 2
params.memory = 6

log.info """\
    C E L L _ L I N E S - N F  P I P E L I N E
    ==========================================
    reads            : ${params.reads}
    reference_dir    : ${params.reference}
    outdir           : ${params.outdir}
    cpus             : ${params.cpus}
    memory           : ${params.memory}
    """
    .stripIndent(true)


/*
 *read_pairs_ch = Channel.fromFilePairs(params.reads)
 *read_pairs_ch.view()
 */


process FASTP_PAIRED {
    tag "$sample_ID"
    cpus "${params.cpus}"
    memory "${params.memory}"
    debug true

    publishDir "${params.outdir}/${sample_ID}/1_fastp", mode: "copy"

    input:
        tuple val(sample_ID), path(reads)

    output:
        tuple val(sample_ID), path("fp_trimmed_*"), emit: fp_trimmed_reads
        tuple val(sample_ID), path("*.json")      , emit: fastp_json
        tuple val(sample_ID), path("*.html")      , emit: fastp_html

    script:
    """
    mkdir -p results/${sample_ID}/1_fastp
    fastp \
    -i ${reads[0]} \
    -I ${reads[1]} \
    -o fp_trimmed_${reads[0]} \
    -O fp_trimmed_${reads[1]} \
    -j ${sample_ID}_fp_trimmed.json \
    -h ${sample_ID}_fp_trimmed.html
    """
}

process DECOMPRESS_FASTQ_GZ {
    tag "$sample_ID"
    cpus "${params.cpus}"
    memory "${params.memory}"
    debug true

    publishDir "${params.outdir}/${sample_ID}/1_fastp", mode: "copy"

    input:
        tuple val(sample_ID), path(fp_trimmed_reads)

    output:
        tuple val(sample_ID), path("d_*"), emit: d_fp_trimmed_reads
    
    script:
    """
    gzip -dc ${fp_trimmed_reads[0]} > d_fp_trimmed_${sample_ID}_1.fq
    gzip -dc ${fp_trimmed_reads[1]} > d_fp_trimmed_${sample_ID}_2.fq
    """
}

process BWA_MEM {
    tag "$sample_ID"
    cpus "${params.cpus}"
    memory "${params.memory}"
    debug true

    publishDir "${params.outdir}/${sample_ID}/2_bwa", mode: "copy"

    input:
        path(reference)
        tuple val(sample_ID), path(d_fp_trimmed_reads)
   
    output:
        tuple val(sample_ID), path("*.bam"), emit: bam
     
    script:
    def idxbase = reference[0].baseName

    """
    mkdir -p results/${sample_ID}/2_bwa
    bwa mem -t ${task.cpus} ${idxbase} ${d_fp_trimmed_reads[0]} ${d_fp_trimmed_reads[1]} | samtools view -Shb -o ./${sample_ID}.bam
    """
}

process SAM_SORT_BAM {
    tag "$sample_ID"
    cpus "${params.cpus}"
    memory "${params.memory}"
    debug true

    publishDir "${params.outdir}/${sample_ID}/2_bwa", mode: "copy"

    input:
        tuple val(sample_ID), path(bam)

    output:
        tuple val(sample_ID), path("*.bam"), emit: sorted_bam

    script:
    """
    samtools sort ${bam} -o sorted_${sample_ID}.bam
    """
}

process SAM_INDEX {
    tag "$sample_ID"
    cpus "${params.cpus}"
    memory "${params. memory}"
    debug true

    publishDir "${params.outdir}/${sample_ID}/2_bwa", mode: "copy"

    input:
        tuple val(sample_ID), path(sorted_bam)

    output:
        tuple val(sample_ID), path("*.bai"), emit: sam_bai

    script:
    """
    samtools index $sorted_bam
    """
}

// Channels
read_pairs_ch = Channel.fromFilePairs(params.reads)
//reference_ch = Channel.fromPath(params.reference)
bwa_index = file('/data/data_AG_Glen/reference_data/Homo_sapiens.GRCh38.dna.primary_assembly.fasta.{,amb,ann,bwt,pac,sa}')

workflow {
    FASTP_PAIRED(read_pairs_ch)
    DECOMPRESS_FASTQ_GZ(FASTP_PAIRED.out.fp_trimmed_reads)
    BWA_MEM(bwa_index,DECOMPRESS_FASTQ_GZ.out.d_fp_trimmed_reads)
    SAM_SORT_BAM(BWA_MEM.out.bam)
    SAM_INDEX(SAM_SORT_BAM.out.sorted_bam)
    SAM_INDEX.out.sam_bai.view()    
}

