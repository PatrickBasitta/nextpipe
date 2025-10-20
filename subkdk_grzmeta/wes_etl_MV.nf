nextflow.enable.dsl=2

process extract_id_filepath {
   
    tag "${target_dir}"
    debug true
    //errorStrategy 'ignore'
    cache 'lenient'

    input:
      path(target_dir)
      val(NextSeq_data_dir)
      val(NovaSeq_data_dir)
      val(nxf_outputdir)

    output:
      path("id_xlsx_paths.csv"), emit: id_xlsx
      path("fastq_paths.csv")  , emit: id_fastqs
      path("bam_path.csv")     , emit: id_bam
      path("vcf_path.csv")     , emit: id_vcf
      path("patient_id.csv")   , emit: id_patient

    script:
    """
    wes_extract_id_filepath.py \\
        --target_dir_mvwes ${target_dir} \\
        --nextseq_data_dir ${NextSeq_data_dir} \\
        --novaseq_data_dir ${NovaSeq_data_dir} \\
        --nxf_outputdir ${nxf_outputdir} \\
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
    //cpus 16
    cache 'lenient'
    errorStrategy { task.exitStatus in 250..253 ? 'terminate' : 'retry' }
    maxRetries 4

    input:
        tuple val(sample_id), path(read1), path(read2)
        val(grz_submission_dir)

    output:
        tuple val(sample_id), path("${read1.getSimpleName()}_${read2.getSimpleName()}_fq_sha256sum.json")        , emit: sha256sum_fqs
        tuple val(sample_id), path("${read1.getSimpleName()}_${read2.getSimpleName()}_fq_bytesize.json")      , emit: bytesize_fqs
        tuple val(sample_id), path("${read1.getSimpleName()}.fastp_1.fq.gz")                                  , emit: fastp_1_fq
        tuple val(sample_id), path("${read2.getSimpleName()}.fastp_2.fq.gz")                                  , emit: fastp_2_fq
        tuple val(sample_id), path("${read1.getSimpleName()}_${read2.getSimpleName()}.fastp.json")            , emit: fastp_json
        tuple val(sample_id), path("${read1.getSimpleName()}_${read2.getSimpleName()}.fastp.html")            , emit: fastp_html
        tuple val(sample_id), path("${read1.getSimpleName()}_${read2.getSimpleName()}.fastp.log")                , emit: fastp_log    

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
        --thread 16 \\
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
    //cpus 8
    //memory "8 GB"

    input:
    tuple val(sample_id), path(bam), path(bai)
    val(grz_submission_dir)
    val(wes_bedfile)
    
    output:
    tuple val(sample_id), path("${bam.getSimpleName()}_depth.stats"), emit: sam_qc
    tuple val(sample_id), path("${bam.getSimpleName()}_bam_qc.csv") , emit: bam_qc
    tuple val(sample_id), path("${bam.getSimpleName()}_bam.json"), emit: bam_json
            
    script:
    def awk1 = "awk -v OFS=',' '{sum+=\$3; ++n; if(\$3>max){max=\$3}; if(\$3<min||min==0){min=\$3};if(\$3>=30){c++}} END{print min, max, sum/n, c/n}'"
    def awk2 = "awk -v OFS=',' '{sum+=\$3; ++n; if(\$3>max){max=\$3}; if(\$3<min||min==0){min=\$3};if(\$3>=100){c++}} END{print min, max, sum/n, c/n}'"
    def awk3 = "awk -v OFS=',' '{num=num?num OFS s1 \$1 s1:s1 \$1 s1} {file=file?file OFS s1 \$2 s1:s1 \$2 s1} END{print num file}'"
    def sample_name = "${bam.getSimpleName()}"
    def  normal_pattern = "N"
    def cmd1 = (sample_name =~ normal_pattern) ? awk1 : awk2
    """
     if [ ! -d ${grz_submission_dir}/${sample_id}/files/ ]
    then
      mkdir -p ${grz_submission_dir}/${sample_id}/files/
    fi

    if [ ! -d ${grz_submission_dir}/${sample_id}/metadata/ ]
    then
      mkdir -p ${grz_submission_dir}/${sample_id}/metadata/
    fi
   
    samtools depth -b ${wes_bedfile} ${bam} --threads 8 -o ${bam.getSimpleName()}_depth.stats
    cat ${bam.getSimpleName()}_depth.stats | ${cmd1} > ${bam.getSimpleName()}_qc_cov.csv
    echo "min_cov,max_cov,mean_cov,targets_above_mincov" > header.csv
    cat  header.csv ${bam.getSimpleName()}_qc_cov.csv | ${awk3} > ${bam.getSimpleName()}_bam_qc.csv
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
        val(wes_bedfile)

    output:
        tuple val(sample_id), path("${vcf.getSimpleName()}_vcf_cal_sha256.txt"), emit: cal_sha256_vcf
        tuple val(sample_id), path("${vcf.getSimpleName()}_vcf_sha256sum.csv") , emit: sha256sum_vcf
        tuple val(sample_id), path("${vcf.getSimpleName()}_vcf_bytesize.csv")  , emit: bytesize_vcf
        tuple val(sample_id), path("${vcf.getSimpleName()}_vcf_sha256sum.json"), emit: json_sha256sum_vcf
        tuple val(sample_id), path("${vcf.getSimpleName()}_vcf_bytesize.json") , emit: json_bytesize_vcf
        tuple path("bed_sha256sum.json"), path("bed_bytesize.json")            , emit: json_sha256_size_bed

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

    sha256sum ${wes_bedfile} > bed_cal_sha256.txt
    cat bed_cal_sha256.txt | ${awk} > bed_sha256sum.csv
    cat bed_sha256sum.csv | jq -Rsn '
                              {"bedfile_checksum":
                                 [inputs
                                  | . / "\n"
                                  | (.[] | select(length > 0) | . / ",") as \$input
                                  | {"file": \$input[1], "fileChecksum": \$input[0]}]}
                                    ' > bed_sha256sum.json

    wc -c ${wes_bedfile} > bed_cal_bytesize.txt
    cat bed_cal_bytesize.txt | ${awk} > bed_bytesize.csv
    cat bed_bytesize.csv | jq -Rsn '
                             {"bed_bytesize":
                                [inputs
                                 | . / "\n"
                                 | (.[] | select(length > 0) | . / ",") as \$input
                                 | {"file":  \$input[1], "fileByteSize": \$input[0]}]}
                                   ' >  bed_bytesize.json

    if [ ! -f ${grz_submission_dir}/${sample_id}/files/${wes_bedfile} ]
    then
      cp ${wes_bedfile} ${grz_submission_dir}/${sample_id}/files/
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
    wes_extract_patient_data.py \\
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

    publishDir(path: "${params.outdir}/${sample_id}/", mode: "copy")

    input:
        tuple val(sample_id), path(xlsx), path(fastp_json_normal), path(sha256sum_fqs_normal), path(bytesize_fqs_normal), path(fastp_json_tumor), path(sha256sum_fqs_tumor), path(bytesize_fqs_tumor), path(bam_json_normal), path(bam_json_tumor), path(json_sha256sum_vcf), path(json_bytesize_vcf), path(patient_data)
        tuple path(sha256sum_bed), path(size_bed)
        val(hgnc)

    output:
        tuple val(sample_id), path("${sample_id}_submit.json"), emit: final_json

    script:
    """
    wes_json_maker.py \\
        --sample_id ${sample_id} \\
        --xlsx_path ${xlsx} \\
        --fastp_json_normal ${fastp_json_normal} \\
        --fq_sha256_json_normal ${sha256sum_fqs_normal} \\
        --fq_bytes_json_normal ${bytesize_fqs_normal} \\
        --fastp_json_tumor ${fastp_json_tumor} \\
        --fq_sha256_json_tumor ${sha256sum_fqs_tumor} \\
        --fq_bytes_json_tumor ${bytesize_fqs_tumor} \\
        --bam_json_normal ${bam_json_normal} \\
        --bam_json_tumor ${bam_json_tumor} \\
        --vcf_sha256_json ${json_sha256sum_vcf} \\
        --vcf_bytes_json ${json_bytesize_vcf} \\
        --patient_data_json ${patient_data} \\
        --bed_sha256_json ${sha256sum_bed} \\
        --bed_bytes_json ${size_bed} \\
        --hgnc ${hgnc}
    """ 
}

workflow wes_ETL_subKDK_grzSubmissionPreparation {

    take:
       target_dir_mvwes_ch
       NextSeq_data_dir_ch
       NovaSeq_data_dir_ch
       nxf_outputdir_ch
       grz_submission_dir_ch
       hgnc_ch
       wes_bedfile_ch       
       
    main:

       extract_id_filepath(target_dir_mvwes_ch,NextSeq_data_dir_ch,NovaSeq_data_dir_ch,nxf_outputdir_ch)
        
       // make samplesheets
       // channel for make_json process
       id_xlsx_ch = extract_id_filepath.out.id_xlsx | splitCsv(header: true) | map { row -> [row.patient_id,
                                                                                             file(row.xlsx_path) ]}
       
       // channel for process_fastp
       id_fastqs_ch = extract_id_filepath.out.id_fastqs | splitCsv(header: true) | map { row ->
                      [[ row.patient_id,
                           file(row.fq_R1_path_n),
                           file(row.fq_R2_path_n),
                           ],
                      [ row.patient_id,
                           file(row.fq_R1_path_t),
                           file(row.fq_R2_path_t),
                           ]
                       ]
                       } | flatMap{ it -> [it[0], it[1]] }

       // id_fastqs_ch.view()
       
       // channel for process_bamfile
       id_bam_ch = extract_id_filepath.out.id_bam | splitCsv(header: true) | map { row -> 
                      [[row.patient_id,
                          file(row.bam_path_n),
                          file(row.bai_path_n),
                          ],
                      [ row.patient_id,
                          file(row.bam_path_t),
                          file(row.bai_path_t),
                          ]
                      ]
                      } | flatMap{ it -> [it[0], it[1]] }

       //id_bam_ch.view()

       // channel for process_vcf
       id_vcf_ch = extract_id_filepath.out.id_vcf | splitCsv(header: true) | map { row -> [row.patient_id,
                                                                                           file(row.vcf_path) ]}
      
       // channel for patient_id
       id_patient_ch = extract_id_filepath.out.id_patient | splitCsv(header: true) | map { row -> [row.patient_id]}
       //id_patient_ch.view()
                                                                                   
       fastq_out = process_fastqs(id_fastqs_ch,grz_submission_dir_ch)
       bam_out = process_bamfile(id_bam_ch,grz_submission_dir_ch,wes_bedfile_ch)
       vcf_bed_out = process_vcf_bedfile(id_vcf_ch,grz_submission_dir_ch,wes_bedfile_ch)

       patient_id = id_patient_ch.flatten()
       patient_data = extract_patient_data(patient_id)

       //fastq_out.fastp_json.view()
       //fastq_out.sha256sum_fqs.view()
       //fastq_out.bytesize_fqs.view()
       joined_fq_data = fastq_out.fastp_json.join(fastq_out.sha256sum_fqs, by:0).join(fastq_out.bytesize_fqs, by:0)
       //joined_fq_data.view()

       //join paired normal tumor and sort N first T last
       paired_fq_data = joined_fq_data.groupTuple(by:0,sort:true).flatten().collect()
       //paired_fq_data.view()
       sorted_normal_tumor_paired_fq_data = paired_fq_data.map{it ->  
                                                                def key = it[0]                        
                                                                def normal_lst = []
                                                                def tumor_lst = []
                                                                // this logic is due that 6 files are expected (3x N, 3x T)
                                                                (it[1].name.findAll {it.contains('N')}) ? normal_lst.add(it[1]) : tumor_lst.add(it[1])
                                                                (it[2].name.findAll {it.contains('N')}) ? normal_lst.add(it[2]) : tumor_lst.add(it[2])
                                                                (it[3].name.findAll {it.contains('N')}) ? normal_lst.add(it[3]) : tumor_lst.add(it[3])
                                                                (it[4].name.findAll {it.contains('N')}) ? normal_lst.add(it[4]) : tumor_lst.add(it[4])
                                                                (it[5].name.findAll {it.contains('N')}) ? normal_lst.add(it[5]) : tumor_lst.add(it[5])
                                                                (it[6].name.findAll {it.contains('N')}) ? normal_lst.add(it[6]) : tumor_lst.add(it[6])
                                                                // output                                                                
                                                                [key,normal_lst,tumor_lst]}.flatten().collect()
                                                             
       //sorted_normal_tumor_paired_fq_data.view()
       
       // join bam json files and sort N first T last       
       bam_data = bam_out.bam_json
       paired_bam_data = bam_data.groupTuple(by:0,sort:true).flatten().collect()
       //paired_bam_data.view()
       
       sorted_normal_tumor_paired_bam_data = paired_bam_data.map{it ->
                                                                def key = it[0]
                                                                def normalbam_lst = []
                                                                def tumorbam_lst = []
                                                                // this logic is due that 1 files are expected (1x N, 1x T)
                                                                (it[1].name.findAll {it.contains('N')}) ? normalbam_lst.add(it[1]) : tumorbam_lst.add(it[1])
                                                                (it[2].name.findAll {it.contains('N')}) ? normalbam_lst.add(it[2]) : tumorbam_lst.add(it[2])
                                                                // output
                                                                [key,normalbam_lst,tumorbam_lst]}.flatten().collect()

       sorted_normal_tumor_paired_bam_data.view()
       
       // join vcf data
       joined_vcf_data = vcf_bed_out.json_sha256sum_vcf.join(vcf_bed_out.json_bytesize_vcf,by:0)
       //joined_vcf_data.view()
       //vcf_bed_out.json_sha256sum_bed.view()
       //paired_vcf_data = joined_vcf_data.groupTuple(by:0,sort:true).flatten().collect()
       //paired_vcf_data.view()

       // bedfile
       wes_bed_ch = vcf_bed_out.json_sha256_size_bed
       // wes_bed_ch.view()
        
       // join data for json
       data_for_json = id_xlsx_ch.join(sorted_normal_tumor_paired_fq_data,by:0).join(sorted_normal_tumor_paired_bam_data,by:0).join(joined_vcf_data,by:0).join(patient_data.meta_data,by:0) 
       data_for_json.view()
       make_json(data_for_json,wes_bed_ch,hgnc_ch)
       

    //emit:
    //    out = json_out.final_json

}

workflow {

    // set channels
    target_dir_mvwes_ch = Channel.fromPath(params.target_dir_mvwes, type: "dir")
    NextSeq_data_dir_ch = Channel.value(params.NextSeq_data_dir)
    NovaSeq_data_dir_ch = Channel.value(params.NovaSeq_data_dir)
    nxf_outputdir_ch = Channel.value(params.nxf_outputdir)
    //pan_bedfile_ch = Channel.value(params.bedfile)
    grz_submission_dir_ch = Channel.value(params.grz_submission_dir)
    hgnc_ch = Channel.value(params.hgnc)
    wes_bedfile_ch = Channel.value(params.wes_bedfile)
    //grz_bed_ch = Channel.value(params.grz_bed)
    //outdir_ch = Channel.value(params.outdir)

    wes_ETL_subKDK_grzSubmissionPreparation(target_dir_mvwes_ch,NextSeq_data_dir_ch,NovaSeq_data_dir_ch,nxf_outputdir_ch,grz_submission_dir_ch,hgnc_ch,wes_bedfile_ch)
}
