nextflow.enable.dsl=2

process extract_id_filepath {
   
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
    
    conda "bioconda::fastp=1.0.1"

    memory = { Math.max(16, (task.attempt * read1.size() * 0.2 / 1000000000).toDouble()) .GB }
    cpus 16
    cache 'lenient'
    errorStrategy { task.exitStatus in 250..253 ? 'terminate' : 'retry' }
    maxRetries 4

    input:
        tuple val(sample_id), path(read1), path(read2)
        val(grz_submission_dir)

    output:
        //tuple val(sample_id), path("${sample_id}_fq_sha256sum.json")        , emit: sha256sum_fqs
        //tuple val(sample_id), path("${sample_id}_fq_bytesize.json")         , emit: bytesize_fqs
        //tuple val(sample_id), path("${read1.getSimpleName()}.fastp_1.fq.gz"), emit: fastp_1_fq
        //tuple val(sample_id), path("${read2.getSimpleName()}.fastp_2.fq.gz"), emit: fastp_2_fq
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
    fastp \\
        --in1 ${read1} \\
        --in2 ${read2} \\
        --out1 ${read1.getSimpleName()}.fastp_1.fq.gz \\
        --out2 ${read2.getSimpleName()}.fastp_2.fq.gz \\
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

process process_bamfile {

    tag "${sample_id}"
    debug true
    //errorStrategy 'ignore'

    conda "bioconda::samtools=1.22.1"
    cache 'lenient'
    cpus 8
    //memory "8 GB"

    input:
    tuple val(sample_id), path(bam)
    val(grz_submission_dir)
    val(pan_bedfile)
    
    output:
    tuple val(sample_id), path("${bam.getSimpleName()}_depth.stats"), emit: sam_qc
    tuple val(sample_id), path("${bam.getSimpleName()}_bam_qc.csv") , emit: bam_qc
    tuple val(sample_id), path("${bam.getSimpleName()}_bam.json")   , emit: bam_json
                
    script:
    def awk1 = "awk -v OFS=',' '{sum+=\$3; ++n; if(\$3>max){max=\$3}; if(\$3<min||min==0){min=\$3};if(\$3>=100){c++}} END{print min, max, sum/n, c/n}'" 
    def awk2 = "awk -v OFS=',' '{num=num?num OFS s1 \$1 s1:s1 \$1 s1} {file=file?file OFS s1 \$2 s1:s1 \$2 s1} END{print num file}'"
    """   
    samtools sort ${bam} -o sorted_${bam}
    samtools index sorted_${bam} sorted_${bam}.bai
    samtools depth -b ${pan_bedfile} sorted_${bam} --threads ${task.cpus} -o ${bam.getSimpleName()}_depth.stats
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
                              ' > ${bam.getSimpleName()}_bam.json
    
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
  
    output:
        tuple val(sample_id), path("${sample_id}_submit.json"), emit: final_json

    script:
    """
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
       bam_out = process_bamfile(id_bam_ch,sub_dir,pan_bedfile_ch)
       vcf_bed_out = process_vcf_bedfile(id_vcf_ch,sub_dir,pan_bedfile_ch)
       
       fq_data = fastq_out.fastp_out
       //fq_data.view()

       bam_data = bam_out.bam_json
       //bam_data.view()

       joined_vcf_data_bed = vcf_bed_out.json_sha256sum_vcf.join(vcf_bed_out.json_bytesize_vcf,by:0).join(vcf_bed_out.json_sha256_size_bed,by:0)
       //joined_vcf_data_bed.view()
   
       // join data for json       
       data_for_json = id_xlsx_ch.join(fq_data,by:0).join(bam_data,by:0).join(joined_vcf_data_bed,by:0).join(patient_data.meta_data,by:0)
       data_for_json.view()
       make_json(data_for_json)
       

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

    pan_ETL_subKDK_grzSubmissionPreparation(target_dir_mvpan_ch,NovaSeq_data_dir_ch,pan_IE_dir_ch,pan_bedfile_ch,grz_submission_dir_ch)
}
