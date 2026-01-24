nextflow.enable.dsl=2

log.info """\
    W G S - C O V E R A G E - P I P E
    ==========================================
    """
    .stripIndent(true)

process FASTP_PAIRED {

    conda "bioconda::fastp=1.0.1"

    tag "$sample_ID"
    cpus 16
    memory = { Math.max(16, (task.attempt * read1.size() * 0.2 / 1000000000).toDouble()) .GB }
    debug true

    publishDir "${params.outdir}/${sample_ID}/1_fastp", mode: "copy"

    input:
        tuple val(sample_ID), path(read1), path(read2)

    output:
        tuple val(sample_ID), path("fp_trimmed_${read1}"), path("fp_trimmed_${read2}"), emit: fp_trimmed_reads
        tuple val(sample_ID), path("${sample_ID}_fp_trimmed.json")      , emit: fastp_json
        tuple val(sample_ID), path("${sample_ID}_fp_trimmed.html")      , emit: fastp_html
        
    script:
    """
    mkdir -p results/${sample_ID}/1_fastp       

    fastp \
    -i ${read1} \
    -I ${read2} \
    -o fp_trimmed_${read1} \
    -O fp_trimmed_${read2} \
    -j ${sample_ID}_fp_trimmed.json \
    -h ${sample_ID}_fp_trimmed.html
    """
}

process DECOMPRESS_FASTQ_GZ {
    tag "$sample_ID"
    cpus 16
    memory = { Math.max(16, (task.attempt * fp_trimmed_read1.size() * 0.2 / 1000000000).toDouble()) .GB }
    debug true

    publishDir "${params.outdir}/${sample_ID}/1_fastp", mode: "copy"

    input:
        tuple val(sample_ID), path(fp_trimmed_read1), path(fp_trimmed_read2)

    output:
        tuple val(sample_ID), path("d_fp_trimmed_${sample_ID}_1.fq"), path("d_fp_trimmed_${sample_ID}_2.fq"),  emit: d_fp_trimmed_reads
    
    script:
    """
    gzip -dc ${fp_trimmed_read1} > d_fp_trimmed_${sample_ID}_1.fq
    gzip -dc ${fp_trimmed_read2} > d_fp_trimmed_${sample_ID}_2.fq
    """
}

process BWA_MEM {

    conda "bioconda::samtools=1.22.1 bioconda::bwa"

    tag "$sample_ID"
    cpus 16
    memory = { Math.max(16, (task.attempt * d_fp_trimmed_read1.size() * 0.2 / 1000000000).toDouble()) .GB }
    debug true

    publishDir "${params.outdir}/${sample_ID}/2_bwa", mode: "copy"

    input:
        path(reference)
        tuple val(sample_ID), path(d_fp_trimmed_read1), path(d_fp_trimmed_read2)
   
    output:
        tuple val(sample_ID), path("${sample_ID}.bam"), emit: bam
     
    script:
    def idxbase = reference[0].baseName

    """
    mkdir -p results/${sample_ID}/2_bwa
    bwa mem -t ${task.cpus} ${idxbase} ${d_fp_trimmed_read1} ${d_fp_trimmed_read2} | samtools view -Shb -o ./${sample_ID}.bam
    """
}

process SAM_SORT_BAM {

    conda "bioconda::samtools=1.22.1"

    tag "$sample_ID"
    cpus 8
    memory = { Math.max(16, (task.attempt * bam.size() * 0.2 / 1000000000).toDouble()) .GB }
    debug true

    publishDir "${params.outdir}/${sample_ID}/2_bwa", mode: "copy"

    input:
        tuple val(sample_ID), path(bam)

    output:
        tuple val(sample_ID), path("sorted_${bam.getSimpleName()}.bam"), emit: sorted_bam

    script:
    """
    samtools sort -@ ${task.cpus} ${bam} -o sorted_${bam.getSimpleName()}.bam
    """
}

process SAM_INDEX_AND_DEPTH {

    conda "bioconda::samtools=1.22.1 bioconda::sambamba bioconda::mosdepth bioconda::picard"

    tag "$sample_ID"
    cpus 8
    memory = { Math.max(16, (task.attempt * sorted_bam.size() * 0.2 / 1000000000).toDouble()) .GB }
    debug true

    publishDir "${params.outdir}/${sample_ID}/3_bamqc", mode: "copy"

    input:
        tuple val(sample_ID), path(sorted_bam)
        val(tmpdir)
        val(grz_bed)

    output:
        tuple val(sample_ID), path("${sorted_bam.getSimpleName()}_markdup_bam_qc.csv") , emit: bam_qc
        tuple val(sample_ID), path("${sorted_bam.getSimpleName()}_markdup_bam.json"), path("${sorted_bam.getSimpleName()}_targeted_markdup_bam.json"), emit: bam_json
        tuple val(sample_ID), path("${sorted_bam.getSimpleName()}_picard_markdup_duplication_metrics.txt"), emit: picard_metrics
        tuple val(sample_ID), path("${sorted_bam.getSimpleName()}_markdup.mosdepth.summary.txt"), emit: mosdepth_summary

    script:
    def awk1 = "awk -v OFS=',' '{sum+=\$4; ++n; if(\$4>=20){c++}} END{print sum/n, c/n}'"
    def awk2 = "awk -v OFS=',' '{sum+=\$4; ++n; if(\$4>=30){c++}} END{print sum/n, c/n}'"
    def awk3 = "awk -v OFS=',' '{num=num?num OFS s1 \$1 s1:s1 \$1 s1} {file=file?file OFS s1 \$2 s1:s1 \$2 s1} END{print num file}'"
    def sample_name = "${sorted_bam.getSimpleName()}"
    def  normal_pattern = "N"
    def cmd1 = (sample_name =~ normal_pattern) ? awk1 : awk2
    """
    samtools index -@ ${task.cpus} ${sorted_bam} -o ${sorted_bam}.bai
    sambamba markdup -t ${task.cpus} --tmpdir ${tmpdir} ${sorted_bam} ${sorted_bam.getSimpleName()}_markdup.bam
    mosdepth --threads ${task.cpus} ${sorted_bam.getSimpleName()}_markdup ${sorted_bam.getSimpleName()}_markdup.bam
    zcat ${sorted_bam.getSimpleName()}_markdup.per-base.bed.gz | ${cmd1} > ${sorted_bam.getSimpleName()}_markdup_qc_cov.csv
    echo "mean_cov,targets_above_mincov" > header.csv
    cat header.csv ${sorted_bam.getSimpleName()}_markdup_qc_cov.csv | ${awk3} > ${sorted_bam.getSimpleName()}_markdup_bam_qc.csv
    qc_info=\$(cat ${sorted_bam.getSimpleName()}_markdup_bam_qc.csv)
    bamfile="${sorted_bam.getSimpleName()}_markdup.bam,"
    full_info=\$bamfile\$qc_info
    echo \$full_info | jq -Rsn '
                         {"bam_qc":
                           [inputs
                            | . / "\n"
                            | (.[] | select(length > 0) | . / ",") as \$input
                            | {"file": \$input[0], "mean_cov": \$input[3], "targets_above_mincov": \$input[4]}]}
                              ' > ${sorted_bam.getSimpleName()}_markdup_bam.json

    mosdepth --by ${grz_bed} --threads ${task.cpus} ${sorted_bam.getSimpleName()}_markdup_regions ${sorted_bam.getSimpleName()}_markdup.bam
    zcat ${sorted_bam.getSimpleName()}_markdup_regions.regions.bed.gz | ${cmd1} > ${sorted_bam.getSimpleName()}_targeted_markdup_qc_cov.csv
    echo "mean_cov,targets_above_mincov" > targeted_header.csv
    cat targeted_header.csv ${sorted_bam.getSimpleName()}_targeted_markdup_qc_cov.csv | ${awk3} > ${sorted_bam.getSimpleName()}_targeted_markdup_bam_qc.csv
    qc_info=\$(cat ${sorted_bam.getSimpleName()}_targeted_markdup_bam_qc.csv)
    bamfile="${sorted_bam.getSimpleName()}_targeted_markdup.bam,"
    full_info=\$bamfile\$qc_info
    echo \$full_info | jq -Rsn '
                         {"targeted_bam_qc":
                           [inputs
                            | . / "\n"
                            | (.[] | select(length > 0) | . / ",") as \$input
                            | {"file": \$input[0], "mean_cov": \$input[3], "targets_above_mincov": \$input[4]}]}
                              ' > ${sorted_bam.getSimpleName()}_targeted_markdup_bam.json

    picard MarkDuplicates \\
        -I ${sorted_bam} \\
        -O ${sorted_bam.getSimpleName()}_picard_markdup.bam \\
        -M ${sorted_bam.getSimpleName()}_picard_markdup_duplication_metrics.txt \\
        --TMP_DIR ${tmpdir}
    """
}

process MOSDEPTH_SUMMARY {

    tag "$sample_ID"
    cpus 16
    memory = { Math.max(16, (task.attempt * sorted_bam.size() * 0.2 / 1000000000).toDouble()) .GB }
    debug true

    publishDir "${params.outdir}/${sample_ID}/3_bamqc", mode: "copy"

    input:
        tuple val(sample_ID), path(mosdepth_summary)

    output:
         tuple val(sample_ID), path("${mosdepth_summary.getSimpleName()}_mosdepth_mean.csv"), emit: mosdepth_result

    script:
    """
    mosdepth_summary.py \\
         --mosdepth_summary ${mosdepth_summary} \\
         --output ${mosdepth_summary.getSimpleName()}_mosdepth_mean.csv
    """

}

// Channels
csv_ch = Channel.fromPath(params.input_csv) | splitCsv(header: true) | map { row-> tuple(row.sampleid, file(row.read1), file(row.read2)) }
reference_ch = Channel.value(params.reference)
bwa_index = file('/projects/reference/sarek_reference/gatk_grch38/Homo_sapiens/GATK/GRCh38/Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64.{alt,amb,ann,bwt,pac,sa}')
grz_bed_ch = Channel.value(params.grz_bedfile)
tmpdir_ch = Channel.value(params.tmpdir)

workflow {
    FASTP_PAIRED(csv_ch)
    DECOMPRESS_FASTQ_GZ(FASTP_PAIRED.out.fp_trimmed_reads)
    BWA_MEM(bwa_index,DECOMPRESS_FASTQ_GZ.out.d_fp_trimmed_reads)
    SAM_SORT_BAM(BWA_MEM.out.bam)
    SAM_INDEX_AND_DEPTH(SAM_SORT_BAM.out.sorted_bam,tmpdir_ch,grz_bed_ch)
    //MOSDEPTH_SUMMARY(SAM_INDEX_AND_DEPTH.out.mosdepth_summary)    
}

