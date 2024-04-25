nextflow.enable.dsl=2

process GENOMIC_AND_META_DATA_TO_JSON {
     
    debug true

    publishDir "${directory}", mode: "copy"

    input:
    val(directory)

    output:
    path("*.json"), emit: json
    env(prefix), emit: prefix   

    script:
    """
    python /data/basitta/nextpipe/mongoDNPM/DNPMdataset_generator_v1.py --path ${directory}
    """
}
/*
process XLSX_TO_ARCHIVE {

    debug true

    input:
    env(prefix)
    val(target_dir)

    script:
    """
    prefix=${prefix}
    target_dir=${target_dir}
    for file in *.xlsx
    do
       if [[ "${prefix[@]}" =~ "${file%%.*}" ]]
       then
          mv "${file}" "${target_dir}"
       fi
    done
    """
}
*/
workflow DNPM_DATA_GENERATOR {

    take:
      args

    main:
      def path = args.path
      def target_dir = args.target_dir
      
      path_ch = Channel.fromPath(path, checkIfExists: true)
      path_ch.view()
      target_dir_ch = Channel.fromPath(target_dir, checkIfExists: true)
      target_dir_ch.view()
      
     GENOMIC_AND_META_DATA_TO_JSON(path_ch)
     prefix_ch = GENOMIC_AND_META_DATA_TO_JSON.out.prefix.flatten()
     prefix_ch.toList().view()
     //XLSX_TO_ARCHIVE(GENOMIC_AND_META_DATA_TO_JSON.out.prefix.collect(),target_dir_ch)
      
}

workflow {
    def args = params
    DNPM_DATA_GENERATOR(args)
}
