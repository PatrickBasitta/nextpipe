nextflow.enable.dsl=2

process extract_id_jsonpath {
   
    tag "${target_dir}"
    debug true
    cache 'lenient'

    input:
      path(target_dir)
      
    output:
      path("id_json_paths.csv"), emit: id_json_path    

    script:
    """
    extract_id_jsonpath.py \\
        --target_dir_jsons ${target_dir} \\
        --id_json_paths id_json_paths.csv
    """
}

process add_bc_uuid_icdo3 {

    tag "${sample_id}"
    debug true
    cache 'lenient'

    secret 'PASSWORD'
    secret 'USERNAME'
    secret 'ca_cert'
    
    publishDir(path: "${params.json_subkdk_dir}/${sample_id}/", pattern: "${sample_id}_patho_meta.json", mode: "copy")
    publishDir(path: "${submission_dir}/${sample_id}/metadata/", pattern: "metadata.json", mode: "copy")
 
    input:
      tuple val(sample_id), path(json_file), val(pid)
      val(genomic_study_subtype)
      val(library_type)
      val(icdo3_xml)
      val(submission_dir)

    output:
      //path("metadata.json"), emit: grz_json
      //path("${json_file.getSimpleName()}_bc.json"), emit: grz_json
      tuple val(sample_id), path("${sample_id}_patho_meta.json"), emit: patho_json
      tuple val(sample_id), path("metadata.json"), emit: grz_json

    script:
    """
    fetch_research_consent.py --patient-id ${pid} --username \$USERNAME --password \$PASSWORD -o ${json_file.getSimpleName()}_bc.json --ca-cert \$ca_cert
    bc_to_meta.py \\
        --sample_id ${sample_id} \\
        --bc_file ${json_file.getSimpleName()}_bc.json \\
        --meta_json ${json_file} \\
        --genomic_study_subtype ${genomic_study_subtype} \\
        --library_type ${library_type} \\
        --patient_id ${pid} \\
        --icdo3_xml ${icdo3_xml}     
    """

}

process grz_validate_encrypt {
    tag "${sample_id}"
    debug true
    cache 'lenient'

    conda "bioconda::grz-cli" 

    input:
        tuple val(sample_id), path(metadata)
        val(submission_dir)
        val(config_file)

    output:
        ("*")

    script:
    """
    echo ${sample_id} > ${sample_id}.log 2>&1 
    cat ${metadata} | jq '.submission.genomicStudySubtype' >> ${sample_id}.log 2>&1
    grz-cli validate --submission-dir ${submission_dir}/${sample_id} --config-file ${config_file} >> ${sample_id}.log 2>&1 
    grz-cli encrypt --submission-dir ${submission_dir}/${sample_id} --config-file ${config_file}  >> ${sample_id}.log 2>&1
    """
}
    
workflow grzSubmissionPreparation {

    take:
       target_dir_json_ch
       genomic_study_subtype_ch
       library_type_ch
       icdo3_xml_ch
       config_file_ch
       submission_dir_ch

    main:
        extract_id_jsonpath(target_dir_json_ch)

       // make samplesheets
       // channel for make_json process
       id_json_ch = extract_id_jsonpath.out.id_json_path | splitCsv(header: true) | map { row -> [row.sample_id,
                                                                                                  file(row.json_path),
                                                                                                  row.pid ]}

       id_json_ch.view()
       genomic_study_subtype_ch.view()
       //submission_dir_ch = Channel.value(params.submission_dir)
       add_bc_uuid_icdo3(id_json_ch,genomic_study_subtype_ch,library_type_ch,icdo3_xml_ch,submission_dir_ch)
       
       id_meta_ch = add_bc_uuid_icdo3.out.grz_json
       grz_validate_encrypt(id_meta_ch,submission_dir_ch,config_file_ch)
}

workflow {

    // set channels
    target_dir_json_ch = Channel.fromPath(params.target_dir_json, type: "dir")
    genomic_study_subtype_ch = Channel.value(params.genomic_study_subtype)
    library_type_ch = Channel.value(params.library_type)
    icdo3_xml_ch = Channel.value(params.icdo3_xml)
    config_file_ch = Channel.value(params.config_file)
    //submission_dir_ch = Channel.fromPath(params.submission_dir, type: "dir")
    submission_dir_ch = Channel.value(params.submission_dir)
    // subworkflow
    grzSubmissionPreparation(target_dir_json_ch,genomic_study_subtype_ch,library_type_ch,icdo3_xml_ch,config_file_ch,submission_dir_ch)
    
}
