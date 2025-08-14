nextflow.enable.dsl=2

process extract_id_filepath {
   
    tag "${target_dir}"
    debug true
    //errorStrategy 'ignore'

    input:
      path(target_dir)
      val(NextSeq_data_dir)
      val(NovaSeq_data_dir)
      val(pan_IE_dir)

    output:
      path("id_xlsx_paths.csv"), emit: id_xlsx
      path("fastq_paths.csv")  , emit: id_fastqs
      path("bam_path.csv")     , emit: id_bam
      path("vcf_path.csv")     , emit: id_vcf

    script:
    """
    extract_id_filepath.py \\
        --target_dir_mvpan ${target_dir} \\
        --nextseq_data_dir ${NextSeq_data_dir} \\
        --novaseq_data_dir ${NovaSeq_data_dir} \\
        --pan_ie_dir ${pan_IE_dir} \\
        --id_xlsx_paths id_xlsx_paths.csv \\
        --id_fastq_paths fastq_paths.csv \\
        --id_bam_path bam_path.csv \\
        --id_vcf_path vcf_path.csv
    """  

}

process process_fastqs {

    tag "${sample_id}"
    debug true
    //errorStrategy 'ignore'

    conda "bioconda::fastp=1.0.0"

    input:
        tuple val(sample_id), path(read1), path(read2)
        val(grz_submission_dir)

    output:
        tuple val(sample_id), path("${sample_id}_fq_sha256sum.json")        , emit: sha256sum_fqs
        tuple val(sample_id), path("${sample_id}_fq_bytesize.json")         , emit: bytesize_fqs
        tuple val(sample_id), path("${read1.getSimpleName()}.fastp_1.fq.gz"), emit: fastp_1_fq
        tuple val(sample_id), path("${read2.getSimpleName()}.fastp_2.fq.gz"), emit: fastp_2_fq
        tuple val(sample_id), path("${sample_id}.merged.fastq.gz"          ), emit: fastp_merged_fq
        tuple val(sample_id), path("${sample_id}.fastp.json")               , emit: fastp_json
        tuple val(sample_id), path("${sample_id}.fastp.html"),                emit: fastp_html
        tuple val(sample_id), path("${sample_id}.fastp.log"),                 emit: fastp_log    

    script:
    def awk = "awk -v OFS=',' '{num=num?num OFS s1 \$1 s1:s1 \$1 s1} {file=file?file OFS s1 \$2 s1:s1 \$2 s1} END{print num,file}'"
    """
    if [ ! -d ${grz_submission_dir}/${sample_id}/files/ ]
    then
      mkdir -p ${grz_submission_dir}/${sample_id}/files/
    fi

    if [ ! -d ${grz_submission_dir}/${sample_id}/metadata/ ]
    then  
      mkdir -p ${grz_submission_dir}/${sample_id}/metadata/
    fi
    
    sha256sum ${read1} ${read2} > ${sample_id}_fq_cal_sha256.txt
    cat ${sample_id}_fq_cal_sha256.txt | ${awk} > ${sample_id}_fq_sha256sum.csv
    cat ${sample_id}_fq_sha256sum.csv | jq -Rsn '
                                           {"fastq_checksums":
                                             [inputs
                                              | . / "\n"
                                              | (.[] | select(length > 0) | . / ",") as \$input
                                              | {"R1": {"file": \$input[2], "fileChecksum": \$input[0]}}, {"R2": {"file": \$input[3], "fileChecksum": \$input[1]}}]}
                                               ' > ${sample_id}_fq_sha256sum.json
    
    wc -c ${read1} > ${sample_id}_fq_cal_bytesize.txt
    wc -c ${read2} >> ${sample_id}_fq_cal_bytesize.txt
    cat ${sample_id}_fq_cal_bytesize.txt | ${awk} > ${sample_id}_fq_bytesize.csv
    cat ${sample_id}_fq_bytesize.csv | jq -Rsn '
                                      {"fastq_bytesizes":
                                         [inputs
                                          | . / "\n"
                                          | (.[] | select(length > 0) | . / ",") as \$input
                                          | {"R1": {"file":  \$input[2], "fileByteSize": \$input[0]}}, {"R2": {"file":  \$input[3],"fileByteSize": \$input[1]}}]}
                                           ' > ${sample_id}_fq_bytesize.json                                                  
    fastp \\
        --in1 ${read1} \\
        --in2 ${read2} \\
        --out1 ${read1.getSimpleName()}.fastp_1.fq.gz \\
        --out2 ${read2.getSimpleName()}.fastp_2.fq.gz \\
        --json ${sample_id}.fastp.json \\
        --html ${sample_id}.fastp.html \\
        -m --merged_out ${sample_id}.merged.fastq.gz \\
        --detect_adapter_for_pe 2> >(tee ${sample_id}.fastp.log >&2)
     
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

process process_bamfile_bedfile {

    tag "${sample_id}"
    debug true
    //errorStrategy 'ignore'

    conda "bioconda::samtools=1.22.1"

    input:
    tuple val(sample_id), path(bam)
    val(bedfile)
    val(grz_submission_dir)
    
    output:
    tuple val(sample_id), path("${bam.getSimpleName()}_depth.stats"), emit: sam_qc
    tuple val(sample_id), path("${bam.getSimpleName()}_bam_qc.csv") , emit: bam_qc
    tuple val(sample_id), path("${sample_id}_bam.json")             , emit: bam_json
    tuple val(sample_id), path("bed_sha256sum.json")                , emit: json_sha256sum_bed
    tuple val(sample_id), path("bed_bytesize.json")                 , emit: json_bytesize_bed
    
            
    script:
    def awk1 = "awk -v OFS=',' '{sum+=\$3; ++n; if(\$3>max){max=\$3}; if(\$3<min||min==0){min=\$3};if(\$3>=100){c++}} END{print min, max, sum/n, c/n}'" 
    def awk2 = "awk -v OFS=',' '{num=num?num OFS s1 \$1 s1:s1 \$1 s1} {file=file?file OFS s1 \$2 s1:s1 \$2 s1} END{print num file}'"
    def awk3 = "awk -v OFS=',' '{num=num?num OFS s1 \$1 s1:s1 \$1 s1} {file=file?file OFS s1 \$2 s1:s1 \$2 s1} END{print num,file}'"
    """
     if [ ! -d ${grz_submission_dir}/${sample_id}/files/ ]
    then
      mkdir -p ${grz_submission_dir}/${sample_id}/files/
    fi

    if [ ! -d ${grz_submission_dir}/${sample_id}/metadata/ ]
    then
      mkdir -p ${grz_submission_dir}/${sample_id}/metadata/
    fi
    
    samtools sort ${bam} -o sorted_${bam}
    samtools index sorted_${bam} sorted_${bam}.bai
    samtools depth -b ${bedfile} -q 30 sorted_${bam} --threads 8 -o ${bam.getSimpleName()}_depth.stats
    cat ${bam.getSimpleName()}_depth.stats | ${awk1} > ${bam.getSimpleName()}_qc_cov.csv
    echo "min_cov,max_cov,mean_cov,targeted_above_mincov" > header.csv
    cat  header.csv ${bam.getSimpleName()}_qc_cov.csv | ${awk2} > ${bam.getSimpleName()}_bam_qc.csv
    qc_info=\$(cat ${bam.getSimpleName()}_bam_qc.csv)
    bamfile="sorted_${bam}"
    full_info=\$bamfile\$qc_info
    echo \$full_info | jq -Rsn '
                         {"bam_qc":
                           [inputs
                            | . / "\n"
                            | (.[] | select(length > 0) | . / ",") as \$input
                            | {"file": \$input[0], "min_cov": \$input[5], "max_cov": \$input[6], "mean_cov": \$input[7], "target_above_mincov": \$input[8]}]}
                              ' > ${sample_id}_bam.json
    
    sha256sum ${bedfile} > bed_cal_sha256.txt
    cat bed_cal_sha256.txt | ${awk3} > bed_sha256sum.csv
    cat bed_sha256sum.csv | jq -Rsn '
                              {"bedfile_checksum":
                                 [inputs
                                  | . / "\n"
                                  | (.[] | select(length > 0) | . / ",") as \$input
                                  | {"file": \$input[1], "fileChecksum": \$input[0]}]}
                                    ' > bed_sha256sum.json

    wc -c ${bedfile} > bed_cal_bytesize.txt
    cat bed_cal_bytesize.txt | ${awk3} > bed_bytesize.csv
    cat bed_bytesize.csv | jq -Rsn '
                             {"bed_bytesize":
                                [inputs
                                 | . / "\n"
                                 | (.[] | select(length > 0) | . / ",") as \$input
                                 | {"file":  \$input[1], "fileByteSize": \$input[0]}]}
                                   ' >  bed_bytesize.json

    if [ ! -f ${grz_submission_dir}/${sample_id}/files/${bedfile} ]
    then
      cp ${bedfile} ${grz_submission_dir}/${sample_id}/files/
    fi
    """
}

process process_vcf {

    tag "${sample_id}"
    debug true
    //errorStrategy 'ignore'

    input:
    tuple val(sample_id), path(vcf)
    val(grz_submission_dir)

    output:
        tuple val(sample_id), path("${sample_id}_vcf_cal_sha256.txt"), emit: cal_ha256_vcf
        tuple val(sample_id), path("${sample_id}_vcf_sha256sum.csv") , emit: sha256sum_vcf
        tuple val(sample_id), path("${sample_id}_vcf_bytesize.csv")  , emit: bytesize_vcf
        tuple val(sample_id), path("${sample_id}_vcf_sha256sum.json"), emit: json_sha256sum_vcf
        tuple val(sample_id), path("${sample_id}_vcf_bytesize.json") , emit: json_bytesize_vcf

    script:
    def awk = "awk -v OFS=',' '{num=num?num OFS s1 \$1 s1:s1 \$1 s1} {file=file?file OFS s1 \$2 s1:s1 \$2 s1} END{print num,file}'"
    """
     if [ ! -d ${grz_submission_dir}/${sample_id}/files/ ]
    then
      mkdir -p ${grz_submission_dir}/${sample_id}/files/
    fi

    if [ ! -d ${grz_submission_dir}/${sample_id}/metadata/ ]
    then
      mkdir -p ${grz_submission_dir}/${sample_id}/metadata/
    fi

    sha256sum ${vcf} > ${sample_id}_vcf_cal_sha256.txt
    cat ${sample_id}_vcf_cal_sha256.txt | ${awk} > ${sample_id}_vcf_sha256sum.csv
    cat ${sample_id}_vcf_sha256sum.csv | jq -Rsn '
                                           {"vcf_checksum":
                                             [inputs
                                              | . / "\n"
                                              | (.[] | select(length > 0) | . / ",") as \$input
                                              | {"file": \$input[1], "fileChecksum": \$input[0]}]}
                                               ' > ${sample_id}_vcf_sha256sum.json

    wc -c ${vcf} > ${sample_id}_vcf_cal_bytesize.txt
    cat ${sample_id}_vcf_cal_bytesize.txt | ${awk} > ${sample_id}_vcf_bytesize.csv
    cat ${sample_id}_vcf_bytesize.csv | jq -Rsn '
                                      {"vcf_bytesize":
                                         [inputs
                                          | . / "\n"
                                          | (.[] | select(length > 0) | . / ",") as \$input
                                          | {"file":  \$input[1], "fileByteSize": \$input[0]}]}
                                           ' > ${sample_id}_vcf_bytesize.json

    if [ ! -f ${grz_submission_dir}/${sample_id}/files/${vcf} ]
    then
      cp ${vcf} ${grz_submission_dir}/${sample_id}/files/
    fi
    """
}

process make_json {

    tag "${sample_id}"
    debug true
    //errorStrategy 'ignore'

    conda "conda-forge::pandas=2.3.1 conda-forge::openpyxl=3.1.5"

    input:
        tuple val(sample_id), path(xlsx), path(fastp_json), path(sha256sum_fqs), path(bytesize_fqs), path(bam_json), path(json_sha256sum_bed), path(json_bytesize_bed), path(json_sha256sum_vcf), path(json_bytesize_vcf)
    
    output:
        tuple val(sample_id), path("${sample_id}_submit.json"), emit: final_json

    script:
    """
    json_maker.py \\
        --sample_id ${sample_id} \\
        --xlsx_path ${xlsx} \\
        --fastp_json ${fastp_json} \\
        --fq_256_json ${sha256sum_fqs} \\
        --fq_bytes_json ${bytesize_fqs} \\
        --bam_json ${bam_json} \\
        --vcf_256_json ${json_sha256sum_vcf} \\
        --vcf_bytes_json ${json_bytesize_vcf} \\
        --bed_256_json ${json_sha256sum_bed} \\
        --bed_bytes_json ${json_bytesize_bed}
    """ 
}

workflow pan_ETL_subKDK_grzSubmissionPreparation {

    take:
       target_dir_mvpan_ch
       NextSeq_data_dir_ch
       NovaSeq_data_dir_ch
       pan_IE_dir_ch
       pan_bedfile_ch
       grz_submission_dir_ch
       
    main:

       extract_id_filepath(target_dir_mvpan_ch,NextSeq_data_dir_ch,NovaSeq_data_dir_ch,pan_IE_dir_ch)
        
       // make samplesheets
       // channel for make_json process
       id_xlsx_ch = extract_id_filepath.out.id_xlsx | splitCsv(header: true) | map { row -> [row.patient_id,
                                                                                             file(row.xlsx_path) ]}
       // channel for process_fastp
       id_fastqs_ch = extract_id_filepath.out.id_fastqs | splitCsv(header: true) | map { row -> [row.patient_id,
                                                                                         file(row.fq_R1_path),
                                                                                         file(row.fq_R2_path) ]}
       // channel for process_bamfile_bedfile
       id_bam_ch = extract_id_filepath.out.id_bam | splitCsv(header: true) | map { row -> [row.patient_id,
                                                                                   file(row.bam_path) ]}
                                                                                   // file(row.bai_path) ]}
       // channel for process_vcf
       id_vcf_ch = extract_id_filepath.out.id_vcf | splitCsv(header: true) | map { row -> [row.patient_id,
                                                                                           file(row.vcf_path) ]}

       fastq_out = process_fastqs(id_fastqs_ch,grz_submission_dir_ch)
       bam_bed_out = process_bamfile_bedfile(id_bam_ch,pan_bedfile_ch,grz_submission_dir_ch)
       vcf_out = process_vcf(id_vcf_ch,grz_submission_dir_ch)

       //fastq_out.fastp_json.view()
       //fastq_out.sha256sum_fqs.view()
       //fastq_out.bytesize_fqs.view()
       joined_fq_data = fastq_out.fastp_json.join(fastq_out.sha256sum_fqs, by:0).join(fastq_out.bytesize_fqs, by:0)
       //joined_fq_data.view()
       joined_bam_bed_data = bam_bed_out.bam_json.join(bam_bed_out.json_sha256sum_bed,by:0).join(bam_bed_out.json_bytesize_bed,by:0)
       //joined_bam_bed_data.view()
       joined_vcf_data = vcf_out.json_sha256sum_vcf.join(vcf_out.json_bytesize_vcf,by:0)
       //joined_vcf_data.view()
       data_for_json = id_xlsx_ch.join(joined_fq_data,by:0).join(joined_bam_bed_data,by:0).join(joined_vcf_data,by:0)
       //data_for_json.view()
       make_json(data_for_json)

}
// set channels
target_dir_mvpan_ch = Channel.fromPath(params.target_dir_mvpan, type: "dir")
NextSeq_data_dir_ch = Channel.value(params.NextSeq_data_dir)
NovaSeq_data_dir_ch = Channel.value(params.NovaSeq_data_dir)
pan_IE_dir_ch = Channel.value(params.pan_IE_dir)
pan_bedfile_ch = Channel.value(params.bedfile)
grz_submission_dir_ch = Channel.value(params.grz_submission_dir)

workflow {
    pan_ETL_subKDK_grzSubmissionPreparation(target_dir_mvpan_ch,NextSeq_data_dir_ch,NovaSeq_data_dir_ch,pan_IE_dir_ch,pan_bedfile_ch,grz_submission_dir_ch)
}
