nextflow.enable.dsl=2

process extract_id_filepath {
   
    tag "${target_dir}"
    debug true
    //errorStrategy 'ignore'
    cache 'lenient'

    input:
      path(target_dir)
      val(NovaSeq_data_dir)
      val(vcf_dir)
      val(fastp_js)
      val(bamqc)
      val(labdata)

    output:
      path("id_xlsx_paths.csv"), emit: id_xlsx
      path("fastq_paths.csv")  , emit: id_fastqs
      path("fastp_path.csv")   , emit: id_fastp_js
      path("bamqc_path.csv")   , emit: id_bamqc
      path("labdata_path.csv") , emit: id_labdata
      path("vcf_path.csv")     , emit: id_vcf
      path("patient_id.csv")   , emit: id_patient

    script:
    """
    wes_extract_id_filepath.py \\
        --target_dir_mvwes ${target_dir} \\
        --novaseq_data_dir ${NovaSeq_data_dir} \\
        --vcf_dir ${vcf_dir} \\
        --fastp_js_dir ${fastp_js} \\
        --bamqc_dir ${bamqc} \\
        --labdata_dir ${labdata} \\
        --id_xlsx_paths id_xlsx_paths.csv \\
        --id_fastq_paths fastq_paths.csv \\
        --id_fastp_path fastp_path.csv \\
        --id_bamqc_path bamqc_path.csv \\
        --id_labdata_path labdata_path.csv \\
        --id_vcf_path vcf_path.csv \\
        --id_patient patient_id.csv
    """  

}

process process_fastqs {

    tag "${sample_id}"
    debug true
    
    //conda "bioconda::fastp=1.0.1"
    
    memory = { Math.max(16, (task.attempt * read1.size() * 0.2 / 1000000000).toDouble()) .GB }
    //cpus 16
    cache 'lenient'
    errorStrategy { task.exitStatus in 250..253 ? 'terminate' : 'retry' }
    maxRetries 4

    input:
        tuple val(sample_id), path(read1), path(read2)
        val(grz_submission_dir)

    output:
        tuple val(sample_id), path("${read1.getSimpleName()}_${read2.getSimpleName()}_fq_sha256sum.json"), path("${read1.getSimpleName()}_${read2.getSimpleName()}_fq_bytesize.json"), emit: fastp_out
   
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
    
    sha256sum ${wes_bedfile} > bed_cal_sha256.txt
    cat bed_cal_sha256.txt | ${awk} > bed_sha256sum.csv
    cat bed_sha256sum.csv | jq -Rsn '
                              {"bedfile_checksum":
                                 [inputs
                                  | . / "\n"
                                  | (.[] | select(length > 0) | . / ",") as \$input
                                  | {"file": \$input[1], "fileChecksum": \$input[0]}]}
                                    ' > ${sample_id}_bed_sha256sum.json

    wc -c ${wes_bedfile} > bed_cal_bytesize.txt
    cat bed_cal_bytesize.txt | ${awk} > bed_bytesize.csv
    cat bed_bytesize.csv | jq -Rsn '
                             {"bed_bytesize":
                                [inputs
                                 | . / "\n"
                                 | (.[] | select(length > 0) | . / ",") as \$input
                                 | {"file":  \$input[1], "fileByteSize": \$input[0]}]}
                                   ' >  ${sample_id}_bed_bytesize.json

    if [ ! -f ${grz_submission_dir}/${sample_id}/files/${wes_bedfile} ]
    then
      cp ${wes_bedfile} ${grz_submission_dir}/${sample_id}/files/
    fi
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
        tuple val(sample_id),
            path(xlsx),
            path(fastp_json_normal),
            path(fastp_json_tumor),
            path(sha256sum_fqs_normal),
            path(bytesize_fqs_normal),
            path(sha256sum_fqs_tumor),
            path(bytesize_fqs_tumor),
            path(bam_json_normal),
            path(bam_json_tumor),
            path(json_sha256sum_vcf),
            path(json_bytesize_vcf),
            path(sha256sum_bed),
            path(size_bed),
            path(labdata)
        val(hgnc)
        val(outdir)

    output:
        tuple val(sample_id), path("${sample_id}_submit.json"), emit: final_json
        tuple val(sample_id), path("${sample_id}.log"), emit: id_log

    script:
    def file_names = [
        "${sample_id}",
        "${xlsx.getSimpleName()}", 
        "${fastp_json_normal.getSimpleName()}",
        "${fastp_json_tumor.getSimpleName()}", 
        "${sha256sum_fqs_normal.getSimpleName()}", 
        "${bytesize_fqs_normal.getSimpleName()}", 
        "${sha256sum_fqs_tumor.getSimpleName()}", 
        "${bytesize_fqs_tumor.getSimpleName()}", 
        "${bam_json_normal.getSimpleName()}", 
        "${bam_json_tumor.getSimpleName()}", 
        "${json_sha256sum_vcf.getSimpleName()}", 
        "${json_bytesize_vcf.getSimpleName()}",  
        "${sha256sum_bed.getSimpleName()}", 
        "${size_bed.getSimpleName()}",
        "${labdata.getSimpleName()}",
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

    wes_json_maker.py \\
        --sample_id ${sample_id} \\
        --xlsx_path ${xlsx} \\
        --fastp_json_normal ${fastp_json_normal} \\
        --fastp_json_tumor ${fastp_json_tumor} \\
        --fq_sha256_json_normal ${sha256sum_fqs_normal} \\
        --fq_bytes_json_normal ${bytesize_fqs_normal} \\
        --fq_sha256_json_tumor ${sha256sum_fqs_tumor} \\
        --fq_bytes_json_tumor ${bytesize_fqs_tumor} \\
        --bam_json_normal ${bam_json_normal} \\
        --bam_json_tumor ${bam_json_tumor} \\
        --vcf_sha256_json ${json_sha256sum_vcf} \\
        --vcf_bytes_json ${json_bytesize_vcf} \\
        --bed_sha256_json ${sha256sum_bed} \\
        --bed_bytes_json ${size_bed} \\
        --labdata_json ${labdata} \\
        --hgnc ${hgnc}
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

workflow wes_ETL_subKDK_grzSubmissionPreparation {

    take:
       target_dir_mvwes_ch
       NovaSeq_data_dir_ch
       grz_submission_dir_ch
       hgnc_ch
       wes_bedfile_ch
       fastp_js_ch
       bamqc_ch
       labdata_ch
       vcf_dir_ch
       outdir_ch

    main:

       extract_id_filepath(target_dir_mvwes_ch,NovaSeq_data_dir_ch,vcf_dir_ch,fastp_js_ch,bamqc_ch,labdata_ch)
        
       // make samplesheets
       // channel for make_json process
       id_xlsx_ch = extract_id_filepath.out.id_xlsx | splitCsv(header: true) | map { row -> [row.patient_id,
                                                                                             file(row.xlsx_path) ]}
       
       // channel for process_fastq
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
       
       // channel for process_vcf
       id_vcf_ch = extract_id_filepath.out.id_vcf | splitCsv(header: true) | map { row -> [row.patient_id,
                                                                                           file(row.vcf_path) ]}

       // channel for fastp_json
       id_fastpjs_ch = extract_id_filepath.out.id_fastp_js | splitCsv(header: true) | map { row -> [row.patient_id,
                                                                                           file(row.fastp_js_path_n),
                                                                                           file(row.fastp_js_path_t) ]}
       
       id_fastpjs_ch.view()

       // channel for bamqc_json
       id_bamqc_ch = extract_id_filepath.out.id_bamqc | splitCsv(header: true) | map { row -> [row.patient_id,
                                                                                           file(row.bamqc_js_path_n),
                                                                                           file(row.bamqc_js_path_t) ]}

       id_bamqc_ch.view()

       // channel for labdata
       id_labdata_ch = extract_id_filepath.out.id_labdata | splitCsv(header: true) | map { row -> [row.patient_id,
                                                                                           file(row.labdata_js_path) ]}
       
       id_labdata_ch.view()

       // channel for patient_id
       id_patient_ch = extract_id_filepath.out.id_patient | splitCsv(header: true) | map { row -> [row.patient_id]}
       id_patient_ch.view()

       patient_id = id_patient_ch.flatten()
       grz_dirs_ch =  grz_dirs(patient_id,grz_submission_dir_ch)
       sub_dir = grz_dirs_ch.done.first()
       //sub_dir.view()
       
                                                                               
       fastq_out = process_fastqs(id_fastqs_ch,sub_dir)
       
       vcf_bed_out = process_vcf_bedfile(id_vcf_ch,sub_dir,wes_bedfile_ch)

       // process fq_data       
       fq_data = fastq_out.fastp_out
       
       //join paired normal tumor and sort N first T last
       paired_fq_data = fq_data.groupTuple() | map {it -> [it[0], it[1..-1].flatten()]}
       //paired_fq_data.view()
       flat_paired_fq_data = paired_fq_data | map { it -> 
                                      def key = it[0]
                                      def val1 = it[1][0]
                                      def val2 = it[1][1]
                                      def val3 = it[1][2]
                                      def val4 = it[1][3]
                                      //def val5 = it[1][4]
                                      //def val6 = it[1][5]
                                      //output
                                      //[key,val1,val2,val3,val4,val5,val6]}
                                      [key,val1,val2,val3,val4]}
       //flat_paired_fq_data.view()
       
       sorted_normal_tumor_paired_fq_data = flat_paired_fq_data.map{it ->  
                                                                def key = it[0]                        
                                                                def normal_lst = []
                                                                def tumor_lst = []
                                                                // this logic is due that 6 files are expected (3x N, 3x T)
                                                                (it[1].name.findAll {it.contains('N')}) ? normal_lst.add(it[1]) : tumor_lst.add(it[1])
                                                                (it[2].name.findAll {it.contains('N')}) ? normal_lst.add(it[2]) : tumor_lst.add(it[2])
                                                                (it[3].name.findAll {it.contains('N')}) ? normal_lst.add(it[3]) : tumor_lst.add(it[3])
                                                                (it[4].name.findAll {it.contains('N')}) ? normal_lst.add(it[4]) : tumor_lst.add(it[4])
                                                                //(it[5].name.findAll {it.contains('N')}) ? normal_lst.add(it[5]) : tumor_lst.add(it[5])
                                                                //(it[6].name.findAll {it.contains('N')}) ? normal_lst.add(it[6]) : tumor_lst.add(it[6])
                                                                // output                                                                
                                                                [key,normal_lst,tumor_lst]} | map { it ->
                                                                                              def key_fq = it[0]
                                                                                              def sha_n = it[1][0]
                                                                                              def byte_n = it[1][1]
                                                                                              //def byte_n = it[1][2]
                                                                                              def sha_t = it[2][0]
                                                                                              def byte_t = it[2][1]
                                                                                              //def byte_t = it[2][2]
                                                                                              //output
                                                                                              //[key_fq,fp_n,sha_n,byte_n,fp_t,sha_t,byte_t]}
                                                                                              [key_fq,sha_n,byte_n,sha_t,byte_t]}
                                                             
       sorted_normal_tumor_paired_fq_data.view()
       
       // join vcf data bedfile
       joined_vcf_data_bed = vcf_bed_out.json_sha256sum_vcf
                                 .join(vcf_bed_out.json_bytesize_vcf,by:0)
                                 .join(vcf_bed_out.json_sha256_size_bed,by:0)
       
       //joined_vcf_data_bed.view()
                    
       // join data for json
       data_for_json = id_xlsx_ch
                          .join(id_fastpjs_ch,by:0)
                          .join(sorted_normal_tumor_paired_fq_data,by:0)
                          .join(id_bamqc_ch,by:0)
                          .join(joined_vcf_data_bed,by:0)
                          .join(id_labdata_ch,by:0)
        
       data_for_json.view()
       data_for_json.count().view()
       //#make_json(data_for_json,hgnc_ch,outdir_ch)
       

    //emit:
    //    out = json_out.final_json

}

workflow {

    // set channels
    target_dir_mvwes_ch = Channel.fromPath(params.target_dir_mvwes, type: "dir")
    // NextSeq_data_dir_ch = Channel.value(params.NextSeq_data_dir)
    NovaSeq_data_dir_ch = Channel.value(params.NovaSeq_data_dir)
    //nxf_outputdir_ch = Channel.value(params.nxf_outputdir)
    //pan_bedfile_ch = Channel.value(params.bedfile)
    grz_submission_dir_ch = Channel.value(params.grz_submission_dir)
    hgnc_ch = Channel.value(params.hgnc)
    wes_bedfile_ch = Channel.value(params.wes_bedfile)
    //grz_bed_ch = Channel.value(params.grz_bed)
    fastp_js_ch = Channel.value(params.fastp_js_dir)
    bamqc_ch = Channel.value(params.bamqc_dir)
    labdata_ch = Channel.value(params.labdata_dir)
    vcf_dir_ch = Channel.value(params.vcf_dir)
    outdir_ch = Channel.value(params.outdir)

    wes_ETL_subKDK_grzSubmissionPreparation(target_dir_mvwes_ch,NovaSeq_data_dir_ch,grz_submission_dir_ch,hgnc_ch,wes_bedfile_ch,fastp_js_ch,bamqc_ch,labdata_ch,vcf_dir_ch,outdir_ch)
}
