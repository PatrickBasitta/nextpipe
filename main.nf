params.reads = "/data/data_AG_Glen/input_data/cl_DE11NGSUKBD135765_5-22Rv1_{1,2}.fq"
params.reference = "/data/data_AG_Glen/reference_data/Homo_sapiens.GRCh38.dna.primary_assembly.fasta"
params.outdir = "results"

log.info """\
    C E L L _ L I N E S - N F   P I P E L I N E
    ===========================================
    reads        : ${params.reads}
    reference    : ${params.reference}
    outdir       : ${params.outdir}
    """
    .stripIndent(true)

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }

    tastp_ch = FASTP_PAIRED(read_pairs_ch)
}

    .set { read_pairs_ch }

/*
 * STEP 1 - FASTP
 */

process FASTP_PAIRED {
    cpus 2
    memory 6
    tag "{name}"

    input:
        tuple name, file(x) from read_pairs_ch

    output:
    path outdir

    script:
    """
    mkdir fastq_trimmed
    fastq \
    -i ${x[0]} -I ${x[1]} \
    -o fastp/trimmed/trim_${x[0]} -O fastq_trimmed/trim_${[1]} }
    -j ${name}_fastq.json
    """     
}

