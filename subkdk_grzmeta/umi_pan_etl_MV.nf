nextflow.enable.dsl=2

process extract_id_filepath {

    cpus 1    
   
    tag "${target_dir}"
    debug true
    //errorStrategy 'ignore'

    input:
      path(target_dir)
      val(NovaSeq_data_dir)
      val(pan_IE_dir)

    output:
      path("id_xlsx_paths.csv"), emit: id_xlsx
      path("fastq_paths.csv")  , emit: id_fastqs
      path("bam_path.csv")     , emit: id_bam
      path("vcf_path.csv")     , emit: id_vcf
      path("patient_id.csv")   , emit: id_patient

    script:
    """
    pan_extract_id_filepath.py \\
        --target_dir_mvpan ${target_dir} \\
        --novaseq_data_dir ${NovaSeq_data_dir} \\
        --pan_ie_dir ${pan_IE_dir} \\
        --id_xlsx_paths id_xlsx_paths.csv \\
        --id_fastq_paths fastq_paths.csv \\
        --id_bam_path bam_path.csv \\
        --id_vcf_path vcf_path.csv \\
        --id_patient patient_id.csv
    """  

}

process process_fastqs {

    tag "${sample_id}"
    debug true
    
    conda "bioconda::fastp=1.0.1 bioconda::umi_tools"

    memory = { Math.max(16, (task.attempt * read1.size() * 0.2 / 1000000000).toDouble()) .GB }
    cpus 8
    cache 'lenient'
    //errorStrategy { task.exitStatus in 250..253 ? 'terminate' : 'retry' }
    //maxRetries 4

    input:
        tuple val(sample_id), path(read1), path(read2)
        val(grz_submission_dir)

    output:
        //tuple val(sample_id), path("${sample_id}_fq_sha256sum.json")        , emit: sha256sum_fqs
        //tuple val(sample_id), path("${sample_id}_fq_bytesize.json")         , emit: bytesize_fqs
        tuple val(sample_id), path("${read1.getSimpleName()}.fastp_1.fq.gz"), path("${read2.getSimpleName()}.fastp_2.fq.gz"), emit: fq_trimmed_reads
        //tuple val(sample_id), path("${sample_id}.merged.fastq.gz"          ), emit: fastp_merged_fq
        //tuple val(sample_id), path("${sample_id}.fastp.json")               , emit: fastp_json
        //tuple val(sample_id), path("${sample_id}.fastp.html"),                emit: fastp_html
        //tuple val(sample_id), path("${sample_id}.fastp.log"),                 emit: fastp_log    
        tuple val(sample_id), path("${read1.getSimpleName()}_${read2.getSimpleName()}.fastp.json"), path("${read1.getSimpleName()}_${read2.getSimpleName()}_fq_sha256sum.json"), path("${read1.getSimpleName()}_${read2.getSimpleName()}_fq_bytesize.json"), emit: fastp_out

    script:
    def awk = "awk -v OFS=',' '{num=num?num OFS s1 \$1 s1:s1 \$1 s1} {file=file?file OFS s1 \$2 s1:s1 \$2 s1} END{print num,file}'"
    """ 
    sha256sum ${read1} ${read2} > ${read1.getSimpleName()}_${read2.getSimpleName()}_fq_cal_sha256.txt
    cat ${read1.getSimpleName()}_${read2.getSimpleName()}_fq_cal_sha256.txt | ${awk} > ${read1.getSimpleName()}_${read2.getSimpleName()}_fq_sha256sum.csv
    cat ${read1.getSimpleName()}_${read2.getSimpleName()}_fq_sha256sum.csv | jq -Rsn '
                                           {"fastq_checksums":
                                             [inputs
                                              | . / "\n"
                                              | (.[] | select(length > 0) | . / ",") as \$input
                                              | {"R1": {"file": \$input[2], "fileChecksum": \$input[0]}}, {"R2": {"file": \$input[3], "fileChecksum": \$input[1]}}]}
                                               ' > ${read1.getSimpleName()}_${read2.getSimpleName()}_fq_sha256sum.json
    
    wc -c ${read1} > ${read1.getSimpleName()}_${read2.getSimpleName()}_fq_cal_bytesize.txt
    wc -c ${read2} >> ${read1.getSimpleName()}_${read2.getSimpleName()}_fq_cal_bytesize.txt
    cat ${read1.getSimpleName()}_${read2.getSimpleName()}_fq_cal_bytesize.txt | ${awk} > ${read1.getSimpleName()}_${read2.getSimpleName()}_fq_bytesize.csv
    cat ${read1.getSimpleName()}_${read2.getSimpleName()}_fq_bytesize.csv | jq -Rsn '
                                      {"fastq_bytesizes":
                                         [inputs
                                          | . / "\n"
                                          | (.[] | select(length > 0) | . / ",") as \$input
                                          | {"R1": {"file":  \$input[2], "fileByteSize": \$input[0]}}, {"R2": {"file": \$input[3],"fileByteSize": \$input[1]}}]}
                                           ' > ${read1.getSimpleName()}_${read2.getSimpleName()}_fq_bytesize.json
    
    umi_tools extract \\
        --bc-pattern2='NNNNNNNNNNNN' \\
        --stdin ${read1} \\
        --stdout ${read1.getSimpleName()}.umi_1.fq.gz \\
        --read2-in ${read2} \\
        --read2-out ${read2.getSimpleName()}.umi_2.fq.gz \\
        -L ${read1.getSimpleName()}_${read2.getSimpleName()}.umi.log
                                                  
    fastp \\
        --in1 ${read1.getSimpleName()}.umi_1.fq.gz \\
        --in2 ${read2.getSimpleName()}.umi_2.fq.gz \\
        --out1 ${read1.getSimpleName()}.fastp_1.fq.gz \\
        --out2 ${read2.getSimpleName()}.fastp_2.fq.gz \\
        --trim_front2 23 \\
        --detect_adapter_for_pe \\
        --trim_poly_g \\
        --cut_mean_quality 20 \\
        --length_required 30 \\
        --thread ${task.cpus} \\
        --json ${read1.getSimpleName()}_${read2.getSimpleName()}.fastp.json \\
        --html ${read1.getSimpleName()}_${read2.getSimpleName()}.fastp.html \\
        2> >(tee ${read1.getSimpleName()}_${read2.getSimpleName()}.fastp.log >&2)
    
    if [ ! -f ${grz_submission_dir}/${sample_id}/files/${read1} ]
    then
      cp ${read1} ${grz_submission_dir}/${sample_id}/files/
    fi
    
    if [ ! -f ${grz_submission_dir}/${sample_id}/files/${read2} ]
    then
      cp ${read2} ${grz_submission_dir}/${sample_id}/files/
    fi

    """
}

process postprocess_fastqs {
    
    tag "${sample_id}"
    cpus 16
    memory = { Math.max(16, (task.attempt * fp_trimmed_read1.size() * 0.2 / 1000000000).toDouble()) .GB }
    debug true

    input:
        tuple val(sample_id), path(fp_trimmed_read1), path(fp_trimmed_read2)
        val(grz_submission_dir)

    output:
        tuple val(sample_id), path("d_fp_trimmed_${fp_trimmed_read1.getSimpleName()}_1.fq"), path("d_fp_trimmed_${fp_trimmed_read2.getSimpleName()}_2.fq"), emit: d_fp_trimmed_reads
    
    script:
    """
    gzip -dc ${fp_trimmed_read1} > d_fp_trimmed_${fp_trimmed_read1.getSimpleName()}_1.fq
    gzip -dc ${fp_trimmed_read2} > d_fp_trimmed_${fp_trimmed_read2.getSimpleName()}_2.fq
    """
}

process bwa_mem {
 
    tag "${sample_id}"
    cpus 16
    memory = { Math.max(16, (task.attempt * fp_trimmed_read1.size() * 0.2 / 1000000000).toDouble()) .GB }
    debug true
    
    conda "bioconda::bwa bioconda::samtools=1.22.1"

    input:
        path(idx)
        tuple val(sample_id), path(fp_trimmed_read1), path(fp_trimmed_read2)
        val(grz_submission_dir)
   
    output:
        tuple val(sample_id), path("${fp_trimmed_read1.getSimpleName()}_sorted.bam"), emit: bwa_mem_bam
     
    script:
    def idxbase = idx[0].baseName
    def lane_id = "lane1"
    def platform_technology = "ILLUMINA"
    def library_prep = "UHS-5000Z-96"
    def platform_unit = "UHS-5000Z-96_L001"
    def sample_name = "${sample_id}"
    def read_group_info = "@RG\\tID:${lane_id}\\tPL:${platform_technology}\\tLB:${library_prep}\\tPU:${platform_unit}\\tSM:${sample_name}"
    """
    bwa mem -t ${task.cpus} -R "${read_group_info}" ${idxbase} ${fp_trimmed_read1} ${fp_trimmed_read2} | samtools view -@ ${task.cpus} -Shb -o ./${fp_trimmed_read1.getSimpleName()}.bam
    samtools sort -@ ${task.cpus} ${fp_trimmed_read1.getSimpleName()}.bam -o ${fp_trimmed_read1.getSimpleName()}_sorted.bam
    """
}

process process_bamfile {

    tag "${sample_id}"
    debug true
    //errorStrategy 'ignore'

    conda "bioconda::samtools=1.22.1 bioconda::sambamba bioconda::fgbio bioconda::umi_tools bioconda::mosdepth"
    cache 'lenient'
    
    cpus { bam.size() > 10.GB ? 8 : 16 }
    memory { bam.size() > 10.GB ? ' 32 GB' : '24 GB' }

    publishDir(path: "${params.outdir}/${sample_id}/coverage", mode: "copy")

    input:
    tuple val(sample_id), path(bam)
    val(grz_submission_dir)
    val(pan_bedfile)
    
    output:
    tuple val(sample_id), path("${bam.getSimpleName()}_depth.stats"), emit: sam_qc
    tuple val(sample_id), path("${bam.getSimpleName()}_bam_qc.csv") , emit: bam_qc
    tuple val(sample_id), path("${bam.getSimpleName().replaceFirst('^fp_trimmed_', '')}_bam.json")   , emit: bam_json
    tuple val(sample_id), path("${bam.getSimpleName()}.bamqc")
                
    script:
    def awk1 = "awk -v OFS=',' '{sum+=\$3; ++n; if(\$3>max){max=\$3}; if(\$3<min||min==0){min=\$3};if(\$3>=100){c++}} END{print min, max, sum/n, c/n}'" 
    def awk2 = "awk -v OFS=',' '{num=num?num OFS s1 \$1 s1:s1 \$1 s1} {file=file?file OFS s1 \$2 s1:s1 \$2 s1} END{print num file}'"
    """
    #samtools view -h ${bam} | samtools view -b -o fixed_${bam}   
    #samtools sort ${bam} -o sorted_${bam}
    samtools index ${bam} ${bam}.bai
    samtools depth -b ${pan_bedfile} ${bam} --threads ${task.cpus} -o ${bam.getSimpleName()}_depth.stats
    cat ${bam.getSimpleName()}_depth.stats | ${awk1} > ${bam.getSimpleName()}_qc_cov.csv
    echo "min_cov,max_cov,mean_cov,targets_above_mincov" > header.csv
    cat  header.csv ${bam.getSimpleName()}_qc_cov.csv | ${awk2} > ${bam.getSimpleName()}_bam_qc.csv
    qc_info=\$(cat ${bam.getSimpleName()}_bam_qc.csv)
    bamfile="${bam},"
    full_info=\$bamfile\$qc_info
    echo \$full_info | jq -Rsn '
                         {"bam_qc":
                           [inputs
                            | . / "\n"
                            | (.[] | select(length > 0) | . / ",") as \$input
                            | {"file": \$input[0], "min_cov": \$input[5], "max_cov": \$input[6], "mean_cov": \$input[7], "targets_above_mincov": \$input[8]}]}
                              ' > ${bam.getSimpleName().replaceFirst('^d_fp_trimmed_', '')}_bam.json

    # markdup with sambamba does not recognize umis
    mkdir tmp_sambamba
    sambamba markdup -t ${task.cpus} --tmpdir tmp_sambamba ${bam} ${bam.getSimpleName()}_markdup.bam
    
    # bam QC with samtools depth without duplicates  
    samtools index -@ ${task.cpus} ${bam.getSimpleName()}_markdup.bam -o ${bam.getSimpleName()}_markdup.bam.bai
    samtools view -b -F 1024 -@ ${task.cpus} ${bam.getSimpleName()}_markdup.bam > filtered_${bam.getSimpleName()}_markdup.bam
    samtools depth -b ${pan_bedfile} -aa -@ ${task.cpus} -s filtered_${bam.getSimpleName()}_markdup.bam > filtered_${bam.getSimpleName()}_markdup.depth
    cat  filtered_${bam.getSimpleName()}_markdup.depth | awk -v OFS=',' -v threshold=100 '{sum+=\$3; ++n; if(\$3>=threshold){c++}} END {print sum/n, c/n}' > filtered_${bam.getSimpleName()}_markdup.cov
    
    # prepare output samtools depth
    echo "file,mean_cov,targets_above_mincov" > header.csv
    sam_results_markdup=\$(cat filtered_${bam.getSimpleName()}_markdup.cov)
    sam_file_markdup=\$(echo "samtools_depth_markdup")
    sam_file_results_markdup=\$(echo \$sam_file_markdup,\$sam_results_markdup)
    echo \$sam_file_results_markdup > sam_file_results_markdup.csv
    cat header.csv sam_file_results_markdup.csv > ${bam.getSimpleName()}.samtools.markdup.depth

    # bam QC with samtools depth including duplicated (old)
    samtools depth -b ${pan_bedfile} -@ ${task.cpus} ${bam.getSimpleName()}_markdup.bam > with_duplicates_${bam.getSimpleName()}_markdup.depth
    cat with_duplicates_${bam.getSimpleName()}_markdup.depth | awk -v OFS=',' -v threshold=100 '{sum+=\$3; ++n; if(\$3>=threshold){c++}} END {print sum/n, c/n}' > with_duplicates_${bam.getSimpleName()}_markdup.cov
    # prepare output samtools depth with duplicates
    sam_results_dup=\$(cat with_duplicates_${bam.getSimpleName()}_markdup.cov)
    sam_file_dup=\$(echo "samtools_depth_w_dup")
    sam_file_results_dup=\$(echo \$sam_file_dup,\$sam_results_dup)
    echo \$sam_file_results_dup > sam_file_dup_results.csv

    # bam QC with mosdepth (for comparison)
    mosdepth --by ${pan_bedfile} ${bam.getSimpleName()}_markdup ${bam.getSimpleName()}_markdup.bam
    zcat  ${bam.getSimpleName()}_markdup.regions.bed.gz | awk '{len = \$3 - \$2; sum += len * \$4; total += len} END {print sum/total}' > ${bam.getSimpleName()}.mosdepth.depth
    zcat  ${bam.getSimpleName()}_markdup.regions.bed.gz | awk -v mincov=100 '{len = \$3 - \$2; total += len; if (\$4 >= mincov) covered += len} END {print (covered/total)}' > ${bam.getSimpleName()}.mosdepth.ontargets
    paste -d',' ${bam.getSimpleName()}.mosdepth.depth ${bam.getSimpleName()}.mosdepth.ontargets > ${bam.getSimpleName()}.mosdepth.region
    cat ${bam.getSimpleName()}_markdup.mosdepth.summary.txt | awk '\$1=="total_region" {print \$4}' > ${bam.getSimpleName()}.mosdepth.summary
    # prepare output mosdepth region
    mosdepth_region_results=\$(cat ${bam.getSimpleName()}.mosdepth.region)
    mosdepth_region_file=\$(echo "mosdepth_region")
    mosdepth_region_file_results=\$(echo \$mosdepth_region_file,\$mosdepth_region_results)
    echo \$mosdepth_region_file_results > mosdepth_region_file_results.csv
    # prepare output mosdepth summary
    mosdepth_summary_results=\$(cat ${bam.getSimpleName()}.mosdepth.summary)
    mosdepth_summary_file=\$(echo "mosdepth_summary")
    mosdepth_summary_file_results=\$(echo \$mosdepth_summary_file,\$mosdepth_summary_results)
    echo \$mosdepth_summary_file_results > mosdepth_summary_file_results.csv

    # join results   
    cat ${bam.getSimpleName()}.samtools.markdup.depth >> ${bam.getSimpleName()}.bamqc
    cat sam_file_dup_results.csv >> ${bam.getSimpleName()}.bamqc
    cat mosdepth_region_file_results.csv >> ${bam.getSimpleName()}.bamqc
    cat mosdepth_summary_file_results.csv >> ${bam.getSimpleName()}.bamqc 

    # umi aware
    # fgbio
    fgbio CopyUmiFromReadName \\
        -i ${bam}  \\
        -o ${bam.getSimpleName()}_umi_tagged.bam \\
        --field-delimiter="_"
   
    # dedup with umi_tools 
    mkdir tmp_dedup
      
    umi_tools dedup \\
        --paired \\
        --extract-umi-method=tag \\
        --umi-tag=RX \\
        --method=directional \\
        --chimeric-pairs=discard \\
        --unpaired-reads=discard \\
        --mapping-quality 30 \\
        -I ${bam.getSimpleName()}_umi_tagged.bam \\
        -S ${bam.getSimpleName()}_dedup.bam \\
        --temp-dir tmp_dedup \\
        --output-stats=${bam.getSimpleName()}_dedup.stats

    # bam QC with samtools depth without duplicates 
    samtools index -@ ${task.cpus} ${bam.getSimpleName()}_dedup.bam -o ${bam.getSimpleName()}_dedup.bam.bai
    samtools depth -b ${pan_bedfile} -aa -@ ${task.cpus} -s ${bam.getSimpleName()}_dedup.bam > filtered_${bam.getSimpleName()}_dedup.depth
    cat filtered_${bam.getSimpleName()}_dedup.depth | awk -v OFS=',' -v threshold=100 '{sum+=\$3; ++n; if(\$3>=threshold){c++}} END {print sum/n, c/n}' > filtered_${bam.getSimpleName()}_dedup.cov

    # prepare output samtools depth
    # echo "file,mean_cov,targets_above_mincov" > bamqc_header.csv
    sam_results_dedup=\$(cat filtered_${bam.getSimpleName()}_dedup.cov)
    sam_file_dedup=\$(echo "samtools_depth_dedup")
    sam_file_results_dedup=\$(echo \$sam_file_dedup,\$sam_results_dedup)
    echo \$sam_file_results_dedup > sam_file_results_dedup.csv
    cat sam_file_results_dedup.csv > ${bam.getSimpleName()}.samtools.dedup.depth

    # bam QC with mosdepth (for comparison)
    mosdepth --by ${pan_bedfile} ${bam.getSimpleName()}_dedup ${bam.getSimpleName()}_dedup.bam
    zcat  ${bam.getSimpleName()}_dedup.regions.bed.gz | awk '{len = \$3 - \$2; sum += len * \$5; total += len} END {print sum/total}' > ${bam.getSimpleName()}.mosdepth.dedup.depth
    zcat  ${bam.getSimpleName()}_dedup.regions.bed.gz | awk -v mincov=100 '{len = \$3 - \$2; total += len; if (\$5 >= mincov) covered += len} END {print (covered/total)}' > ${bam.getSimpleName()}.mosdepth.dedup.ontargets
    paste -d',' ${bam.getSimpleName()}.mosdepth.dedup.depth ${bam.getSimpleName()}.mosdepth.dedup.ontargets > ${bam.getSimpleName()}.mosdepth.dedup.region
    cat ${bam.getSimpleName()}_dedup.mosdepth.summary.txt | awk '\$1=="total_region" {print \$4}' > ${bam.getSimpleName()}.mosdepth.dedup.summary

    # prepare output mosdepth region
    mosdepth_region_results_dedup=\$(cat ${bam.getSimpleName()}.mosdepth.dedup.region)
    mosdepth_region_file_dedup=\$(echo "mosdepth_region_dedup")
    mosdepth_region_file_results_dedup=\$(echo \$mosdepth_region_file_dedup,\$mosdepth_region_results_dedup)
    echo \$mosdepth_region_file_results_dedup > mosdepth_region_file_results_dedup.csv
    # prepare output mosdepth summary
    mosdepth_summary_results_dedup=\$(cat ${bam.getSimpleName()}.mosdepth.dedup.summary)
    mosdepth_summary_file_dedup=\$(echo "mosdepth_summary_dedup")
    mosdepth_summary_file_results_dedup=\$(echo \$mosdepth_summary_file_dedup,\$mosdepth_summary_results_dedup)
    echo \$mosdepth_summary_file_results_dedup > mosdepth_summary_file_results_dedup.csv
    # join results   
    cat ${bam.getSimpleName()}.samtools.dedup.depth >> ${bam.getSimpleName()}.bamqc
    cat mosdepth_region_file_results_dedup.csv >> ${bam.getSimpleName()}.bamqc
    cat mosdepth_summary_file_results_dedup.csv >> ${bam.getSimpleName()}.bamqc     
    """
}

process process_vcf_bedfile {

    tag "${sample_id}"
    debug true
    //errorStrategy 'ignore'
    cache 'lenient' 

    input:
        tuple val(sample_id), path(vcf)
        val(grz_submission_dir)
        val(pan_bedfile)

    output:
        tuple val(sample_id), path("${sample_id}_vcf_cal_sha256.txt"), emit: cal_ha256_vcf
        tuple val(sample_id), path("${sample_id}_vcf_sha256sum.csv") , emit: sha256sum_vcf
        tuple val(sample_id), path("${sample_id}_vcf_bytesize.csv")  , emit: bytesize_vcf
        tuple val(sample_id), path("${sample_id}_vcf_sha256sum.json"), emit: json_sha256sum_vcf
        tuple val(sample_id), path("${sample_id}_vcf_bytesize.json") , emit: json_bytesize_vcf
        tuple val(sample_id), path("${sample_id}_bed_sha256sum.json"), path("${sample_id}_bed_bytesize.json") , emit: json_sha256_size_bed

    script:
    def awk = "awk -v OFS=',' '{num=num?num OFS s1 \$1 s1:s1 \$1 s1} {file=file?file OFS s1 \$2 s1:s1 \$2 s1} END{print num,file}'"
    """
    sha256sum ${vcf} > ${vcf.getSimpleName()}_vcf_cal_sha256.txt
    cat ${vcf.getSimpleName()}_vcf_cal_sha256.txt | ${awk} > ${vcf.getSimpleName()}_vcf_sha256sum.csv
    cat ${vcf.getSimpleName()}_vcf_sha256sum.csv | jq -Rsn '
                                           {"vcf_checksum":
                                             [inputs
                                              | . / "\n"
                                              | (.[] | select(length > 0) | . / ",") as \$input
                                              | {"file": \$input[1], "fileChecksum": \$input[0]}]}
                                               ' > ${vcf.getSimpleName()}_vcf_sha256sum.json

    wc -c ${vcf} > ${vcf.getSimpleName()}_vcf_cal_bytesize.txt
    cat ${vcf.getSimpleName()}_vcf_cal_bytesize.txt | ${awk} > ${vcf.getSimpleName()}_vcf_bytesize.csv
    cat ${vcf.getSimpleName()}_vcf_bytesize.csv | jq -Rsn '
                                      {"vcf_bytesize":
                                         [inputs
                                          | . / "\n"
                                          | (.[] | select(length > 0) | . / ",") as \$input
                                          | {"file":  \$input[1], "fileByteSize": \$input[0]}]}
                                           ' > ${vcf.getSimpleName()}_vcf_bytesize.json

    if [ ! -f ${grz_submission_dir}/${sample_id}/files/${vcf} ]
    then
      cp ${vcf} ${grz_submission_dir}/${sample_id}/files/
    fi

    sha256sum ${pan_bedfile} > bed_cal_sha256.txt
    cat bed_cal_sha256.txt | ${awk} > bed_sha256sum.csv
    cat bed_sha256sum.csv | jq -Rsn '
                              {"bedfile_checksum":
                                 [inputs
                                  | . / "\n"
                                  | (.[] | select(length > 0) | . / ",") as \$input
                                  | {"file": \$input[1], "fileChecksum": \$input[0]}]}
                                    ' > ${sample_id}_bed_sha256sum.json

    wc -c ${pan_bedfile} > bed_cal_bytesize.txt
    cat bed_cal_bytesize.txt | ${awk} > bed_bytesize.csv
    cat bed_bytesize.csv | jq -Rsn '
                             {"bed_bytesize":
                                [inputs
                                 | . / "\n"
                                 | (.[] | select(length > 0) | . / ",") as \$input
                                 | {"file":  \$input[1], "fileByteSize": \$input[0]}]}
                                   ' >  ${sample_id}_bed_bytesize.json

    if [ ! -f ${grz_submission_dir}/${sample_id}/files/${pan_bedfile} ]
    then
      cp ${pan_bedfile} ${grz_submission_dir}/${sample_id}/files/
    fi
    """
}

process extract_patient_data {

    tag "${sample_id}"
    debug true
    //errorStrategy 'ignore'
    cache 'lenient'

    input:
        val(sample_id)

    output:
        tuple val(sample_id), path("${sample_id}_meta_data.json"), emit: meta_data

    script:
    """
    pan_extract_patient_data.py \\
         --patient_id ${sample_id} \\
         --patient_meta_data ${sample_id}_meta_data.json
    """
}

process make_json {

    tag "${sample_id}"
    debug true
    //errorStrategy 'ignore'
    cache 'lenient'

    conda "conda-forge::pandas=2.3.1 conda-forge::openpyxl=3.1.5 conda-forge::requests"

    publishDir(path: "${outdir}/${sample_id}/", mode: "copy")

    input:
        tuple val(sample_id), path(xlsx), path(fastp_json), path(sha256sum_fqs), path(bytesize_fqs), path(bam_json), path(json_sha256sum_vcf), path(json_bytesize_vcf), path(json_sha256sum_bed), path(json_bytesize_bed), path(patient_data)
        val(outdir)
  
    output:
        tuple val(sample_id), path("${sample_id}_submit.json"), emit: final_json
        tuple val(sample_id), path("${sample_id}.log"), emit: id_log

    script:
    def file_names = [
        "${sample_id}",
        "${xlsx.getSimpleName()}", 
        "${fastp_json.getSimpleName()}", 
        "${sha256sum_fqs.getSimpleName()}", 
        "${bytesize_fqs.getSimpleName()}", 
        "${bam_json.getSimpleName()}", 
        "${json_sha256sum_vcf.getSimpleName()}", 
        "${json_bytesize_vcf.getSimpleName()}",
        "${json_sha256sum_bed.getSimpleName()}",
        "${json_bytesize_bed.getSimpleName()}", 
        "${patient_data.getSimpleName()}" 
        ]
    """
    id=${file_names[0]}
  
    if [[ "\${id}" ==  "${sample_id}" ]]; then
        for file_name in ${file_names.join(" ")}; do
            extracted_id=\$(echo \$file_name | cut -d'_' -f1 | cut -d'-' -f1-2)
            if [[ "\$extracted_id" != "${sample_id}" ]]; then
                echo "sample_id \${id} does not match file \$file_name" >> ${sample_id}.log
                exit 1
            else
                echo "sample_id \${id} does match file \$file_name" >> ${sample_id}.log
            fi
        done
    fi  

    pan_json_maker.py \\
        --sample_id ${sample_id} \\
        --xlsx_path ${xlsx} \\
        --fastp_json ${fastp_json} \\
        --fq_256_json ${sha256sum_fqs} \\
        --fq_bytes_json ${bytesize_fqs} \\
        --bam_json ${bam_json} \\
        --vcf_256_json ${json_sha256sum_vcf} \\
        --vcf_bytes_json ${json_bytesize_vcf} \\
        --bed_256_json ${json_sha256sum_bed} \\
        --bed_bytes_json ${json_bytesize_bed} \\
        --patient_data_json ${patient_data}
    """ 
}

process grz_dirs {

    tag "${sample_id}"
    debug true
    //errorStrategy 'ignore'
    cache 'lenient'

    input:
        val(sample_id)
        val(grz_submission_dir)

    output:
        val(grz_submission_dir), emit: done
    
    script:
    """
    if [ ! -d ${grz_submission_dir}/${sample_id}/files/ ]
    then
      mkdir -p ${grz_submission_dir}/${sample_id}/files/
    fi

    if [ ! -d ${grz_submission_dir}/${sample_id}/metadata/ ]
    then
      mkdir -p ${grz_submission_dir}/${sample_id}/metadata/
    fi
    """
}

workflow pan_ETL_subKDK_grzSubmissionPreparation {

    take:
       target_dir_mvpan_ch
       NovaSeq_data_dir_ch
       pan_IE_dir_ch
       pan_bedfile_ch
       grz_submission_dir_ch
       outdir_ch
       bwa_index_ch
       gatk_GRCh38_ref_ch
       
       
    main:

       extract_id_filepath(target_dir_mvpan_ch,NovaSeq_data_dir_ch,pan_IE_dir_ch)
        
       // make samplesheets
       // channel for make_json process
       id_xlsx_ch = extract_id_filepath.out.id_xlsx | splitCsv(header: true) | map { row -> [row.patient_id,
                                                                                             file(row.xlsx_path) ]}
       // channel for process_fastp
       id_fastqs_ch = extract_id_filepath.out.id_fastqs | splitCsv(header: true) | map { row -> [row.patient_id,
                                                                                         file(row.fq_R1_path),
                                                                                         file(row.fq_R2_path) ]}
       id_fastqs_ch.view()
       // channel for process_bamfile_bedfile
       id_bam_ch = extract_id_filepath.out.id_bam | splitCsv(header: true) | map { row -> [row.patient_id,
                                                                                   file(row.bam_path) ]}
                                                                                   // file(row.bai_path) ]}
       // channel for process_vcf
       id_vcf_ch = extract_id_filepath.out.id_vcf | splitCsv(header: true) | map { row -> [row.patient_id,
                                                                                           file(row.vcf_path) ]}
      
       // channel for patient_id
       id_patient_ch = extract_id_filepath.out.id_patient | splitCsv(header: true) | map { row -> [row.patient_id]}
       id_patient_ch.view()
       
       patient_id = id_patient_ch.flatten()
       grz_dirs_ch =  grz_dirs(patient_id,grz_submission_dir_ch)
       sub_dir = grz_dirs_ch.done.first()
       //sub_dir.view()
       patient_data = extract_patient_data(patient_id)
                                                                                
       fastq_out = process_fastqs(id_fastqs_ch,sub_dir)
       
       // prepare for bwa-mem
       trimmed_fq_reads = fastq_out.fq_trimmed_reads
       //trimmed_fq_reads.view()
       
       //#bwa_mem_ready_fqs = postprocess_fastqs(trimmed_fq_reads,sub_dir)
       bam_new = bwa_mem(bwa_index_ch,trimmed_fq_reads,sub_dir)
       bam_out = process_bamfile(bam_new.bwa_mem_bam,sub_dir,pan_bedfile_ch)
       vcf_bed_out = process_vcf_bedfile(id_vcf_ch,sub_dir,pan_bedfile_ch)
       
       fq_data = fastq_out.fastp_out
       //fq_data.view()

       bam_data = bam_out.bam_json
       //bam_data.view()

       joined_vcf_data_bed = vcf_bed_out.json_sha256sum_vcf.join(vcf_bed_out.json_bytesize_vcf,by:0).join(vcf_bed_out.json_sha256_size_bed,by:0)
       //joined_vcf_data_bed.view()
   
       // join data for json       
       data_for_json = id_xlsx_ch.join(fq_data,by:0).join(bam_data,by:0).join(joined_vcf_data_bed,by:0).join(patient_data.meta_data,by:0)
       //data_for_json.view()
       make_json(data_for_json,outdir_ch)
       

    //emit:
    //    out = json_out.final_json

}

workflow {

    // set channels
    target_dir_mvpan_ch = Channel.fromPath(params.target_dir_mvpan, type: "dir")
    // NextSeq_data_dir_ch = Channel.value(params.NextSeq_data_dir)
    NovaSeq_data_dir_ch = Channel.value(params.NovaSeq_data_dir)
    pan_IE_dir_ch = Channel.value(params.pan_IE_dir)
    pan_bedfile_ch = Channel.value(params.bedfile)
    grz_submission_dir_ch = Channel.value(params.grz_submission_dir)
    outdir_ch = Channel.value(params.outdir)
    //bwa_index_ch = file('/projects/reference/sarek_reference/gatk_grch38/Homo_sapiens/GATK/GRCh38/Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64.{alt,amb,ann,bwt,pac,sa}')
    bwa_index_ch = file(params.bwa_index)
    gatk_GRCh38_ref_ch = Channel.value(params.gatk_GRCh38_ref) 

    pan_ETL_subKDK_grzSubmissionPreparation(target_dir_mvpan_ch,NovaSeq_data_dir_ch,pan_IE_dir_ch,pan_bedfile_ch,grz_submission_dir_ch,outdir_ch,bwa_index_ch,gatk_GRCh38_ref_ch)
}
