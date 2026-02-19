nextflow.enable.dsl=2

process extract_id_bam_filepath {
   
    tag "${target_dir}"
    debug true
    //errorStrategy 'ignore'
    cache 'lenient'
   
    //conda "conda-forge::pandas=2.2.3"
    conda "conda-forge::python=3.9.15 conda-forge::pandas=2.0.3"
    container "biocontainers/mulled-v2-0594c09780adaaa41fe60b1869ba41c8905a0c98:24a8102d6795963b77f04bb83cc82c081e4a2adc-0"

    input:
      path(target_dir)
      
    output:
      path("bam_path.csv")     , emit: id_bam
      path("purple_path.csv")  , emit: id_purple
      
    script:
    """
    id_bam_path.py \\
        --target_dir ${target_dir} \\
        --id_bam_path bam_path.csv \\
        --id_purple_path purple_path.csv
    """  
}

process bam_sort_index {

    tag "${bam.getSimpleName()}"
    //cpus 16
    memory = { Math.max(16, (task.attempt * bam.size() * 0.2 / 1000000000).toDouble()) .GB }
    debug true

    conda "bioconda::samtools=1.22.1"

    //publishDir "${params.outdir}/${sample_ID}", mode: "copy"

    input:
        tuple val(sample_ID), path(bam), path(bai)

    output:
        tuple val(sample_ID), path("sorted_${bam.getSimpleName()}.bam"), path("sorted_${bam.getSimpleName()}.bam.bai"), emit: si_bam

    script:
    """
    samtools sort -@ 8 ${bam} -o sorted_${bam.getSimpleName()}.bam
    samtools index -@ 8 sorted_${bam.getSimpleName()}.bam -o sorted_${bam.getSimpleName()}.bam.bai
    """
}
    
process bam_qc {
    
    tag "${sorted_bam.getSimpleName()}"
    //cpus 16
    memory = { Math.max(16, (task.attempt * sorted_bam.size() * 0.2 / 1000000000).toDouble()) .GB }
    debug true

    conda "bioconda::samtools=1.22.1 bioconda::sambamba bioconda::mosdepth"

    publishDir "${params.outdir}/wes_${sample_ID}/coverage", mode: "copy"

    input:
        tuple val(sample_ID), path(sorted_bam), path(sorted_bai)
        val(wes_bedfile)
        val(tmpdir)
        val(library_type)

    output:
        tuple val(sample_ID), path("${sorted_bam.getSimpleName()}.samtools.depth"), path("${sorted_bam.getSimpleName()}.mosdepth.region"), path("${sorted_bam.getSimpleName()}.mosdepth.summary"), path("${sorted_bam.getSimpleName()}.bamqc")
        tuple val(sample_ID), path("filtered_${sorted_bam.getSimpleName()}_markdup.bam"), emit: filtered_bam

    script:
    def sample_name = "${sorted_bam.getSimpleName()}"
    //def  normal_pattern = "NB"
    //def cmd1 = (sample_name =~ /[normal_pattern]/) ? awk1 : awk2
    //def cmd2 = (sample_name =~ /[normal_pattern]/) ? awk3 : awk4
    """
    # markdup with sambamba
    sambamba markdup -t 8 --tmpdir ${tmpdir} ${sorted_bam} ${sorted_bam.getSimpleName()}_markdup.bam
    # bam QC with samtools
    if [[ "${sample_name}" =~ [NB] && "${library_type}" == "wes" ]]; then
        cov_value=30
    elif [[ "${sample_name}" =~ [T] && "${library_type}" == "wes" ]]; then
        cov_value=100
    elif [[ "${sample_name}" =~ [NB] && "${library_type}" == "wgs" ]]; then
        cov_value=20
    elif [[ "${sample_name}" =~ [T] && "${library_type}" == "wgs" ]]; then
        cov_value=30
    else
        echo "No matching filename or wrong library_type"
        exit 1
    fi
     
    samtools index -@20 ${sorted_bam.getSimpleName()}_markdup.bam -o ${sorted_bam.getSimpleName()}_markdup.bam.bai
    samtools view -b -F 1024 -@ 50 ${sorted_bam.getSimpleName()}_markdup.bam > filtered_${sorted_bam.getSimpleName()}_markdup.bam
    samtools depth -b ${wes_bedfile} -aa -@ 50 -s filtered_${sorted_bam.getSimpleName()}_markdup.bam > filtered_${sorted_bam.getSimpleName()}_markdup.depth 
    cat  filtered_${sorted_bam.getSimpleName()}_markdup.depth | awk -v OFS=',' -v threshold=\$cov_value '{sum+=\$3; ++n; if(\$3>=threshold){c++}} END {print sum/n, c/n}' > filtered_${sorted_bam.getSimpleName()}_markdup.cov
    # prepare output samtools depth
    echo "tool,mean_cov,targets_above_mincov" > header.csv
    sam_results=\$(cat filtered_${sorted_bam.getSimpleName()}_markdup.cov)
    sam_file=\$(echo "samtools_depth")
    sam_file_results=\$(echo \$sam_file,\$sam_results)
    echo \$sam_file_results > sam_file_results.csv
    cat header.csv sam_file_results.csv > ${sorted_bam.getSimpleName()}.samtools.depth
    # bam QC with mosdepth (for comparison)
    mosdepth --by ${wes_bedfile} ${sorted_bam.getSimpleName()}_markdup ${sorted_bam.getSimpleName()}_markdup.bam
    zcat  ${sorted_bam.getSimpleName()}_markdup.regions.bed.gz | awk '{len = \$3 - \$2; sum += len * \$4; total += len} END {print sum/total}' > ${sorted_bam.getSimpleName()}.mosdepth.depth
    zcat  ${sorted_bam.getSimpleName()}_markdup.regions.bed.gz | awk -v mincov=\$cov_value '{len = \$3 - \$2; total += len; if (\$4 >= mincov) covered += len} END {print (covered/total)}' > ${sorted_bam.getSimpleName()}.mosdepth.ontargets
    paste -d',' ${sorted_bam.getSimpleName()}.mosdepth.depth ${sorted_bam.getSimpleName()}.mosdepth.ontargets > ${sorted_bam.getSimpleName()}.mosdepth.region
    cat ${sorted_bam.getSimpleName()}_markdup.mosdepth.summary.txt | awk '\$1=="total_region" {print \$4}' > ${sorted_bam.getSimpleName()}.mosdepth.summary
    # prepare output mosdepth region
    mosdepth_region_results=\$(cat ${sorted_bam.getSimpleName()}.mosdepth.region)
    mosdepth_region_file=\$(echo "mosdepth_region")
    mosdepth_region_file_results=\$(echo \$mosdepth_region_file,\$mosdepth_region_results)
    echo \$mosdepth_region_file_results > mosdepth_region_file_results.csv
    # prepare output mosdepth summary
    mosdepth_summary_results=\$(cat ${sorted_bam.getSimpleName()}.mosdepth.summary)
    mosdepth_summary_file=\$(echo "mosdepth_summary")
    mosdepth_summary_file_results=\$(echo \$mosdepth_summary_file,\$mosdepth_summary_results)
    echo \$mosdepth_summary_file_results > mosdepth_summary_file_results.csv
    # join results   
    cat ${sorted_bam.getSimpleName()}.samtools.depth  >> ${sorted_bam.getSimpleName()}.bamqc
    cat mosdepth_region_file_results.csv >> ${sorted_bam.getSimpleName()}.bamqc
    cat mosdepth_summary_file_results.csv >> ${sorted_bam.getSimpleName()}.bamqc 
    """
}

process index_filtered_bams {

    tag "${sample_ID}"
    cpus 8
    memory = { Math.max(16, (task.attempt * filtered_tumorbam.size() * 0.2 / 1000000000).toDouble()) .GB }
    debug true

    conda "bioconda::samtools=1.22.1"

    input:
        tuple val(sample_ID), path(filtered_normalbam), path(filtered_tumorbam)

    output:
        tuple val(sample_ID), path("${filtered_normalbam}"), path("${filtered_normalbam}.bai"), path("${filtered_tumorbam}"), path("${filtered_tumorbam}.bai"), emit: idx_filtered_bam

    script:
    """
    samtools index -@ 24 ${filtered_normalbam} -o ${filtered_normalbam}.bai
    samtools index -@ 24 ${filtered_tumorbam} -o ${filtered_tumorbam}.bai
    """ 

}

process msisensor_pro {

    tag "${sample_ID}"
    //cpus 16
    memory = { Math.max(16, (task.attempt * idx_filtered_tumorbam.size() * 0.2 / 1000000000).toDouble()) .GB }
    debug true

    conda "bioconda::msisensor-pro=1.2.0"

    publishDir "${params.outdir}/wes_${sample_ID}/msisensor_pro", mode: "copy"

    input:
        tuple val(sample_ID), path(idx_filtered_normalbam), path(idx_filtered_normalbai), path(idx_filtered_tumorbam), path(idx_filtered_tumorbai) 
        val(ref_genome_fasta)
        //val(ref_genome_fai)
        //val(ref_genome_dict)

    output:
        tuple val(sample_ID), path("${idx_filtered_normalbam.getSimpleName()}_${idx_filtered_tumorbam.getSimpleName()}_msi.csv")
        //path("${ref_genome_fasta}_msi.list")
    script:
    """
    msisensor-pro scan \\
        -d ${ref_genome_fasta} \\
        -o msi.list

    msisensor-pro msi \\
        -d msi.list \\
        -n ${idx_filtered_normalbam} \\
        -t ${idx_filtered_tumorbam} \\
        -o ${idx_filtered_normalbam.getSimpleName()}_${idx_filtered_tumorbam.getSimpleName()}_msi.csv
    """
}

workflow postprocessing_wes_custom_oa {

    take:
       target_dir_ch
       wes_bedfile_ch
       tmpdir_ch
       ref_genome_fasta_ch
       ref_genome_fai_ch
       ref_genome_dict_ch
       library_type_ch
       outdir_ch 

    main:

       extract_id_bam_filepath(target_dir_ch)

       // channel for process_bamfile
       id_bam_qc_ch = extract_id_bam_filepath.out.id_bam | splitCsv(header: true) | map { row -> 
                       [[row.patient_id,
                           file(row.normal_bam),
                           file(row.normal_bai),
                           ],
                       [row.patient_id,
                           file(row.tumor_bam),
                           file(row.tumor_bai),
                           ]
                       ]
                       } | flatMap{ it -> [it[0], it[1]] }

       id_bam_qc_ch.view()

       // channel purple
       //id_purple_chord_ch = extract_id_bam_filepath.out.id_purple | splitCsv(header: true) | map { row -> [row.patient_id,file(row.purple_snv),file(row.purple_sv)]}
     
       //id_purple_chord_ch.view()

       // BAM QC
       sorted_and_indexed_bam = bam_sort_index(id_bam_qc_ch)
       bam_qc(sorted_and_indexed_bam.si_bam,wes_bedfile_ch,tmpdir_ch,library_type_ch)
       //bam_qc.out.filtered_bam.view() 
       // join and sort tumor normal bam

       // Prepare fpr msisensor pro
       joined_bams = bam_qc.out.filtered_bam.groupTuple() 
       //joined_bams.view()
       sorted_joined_bams = joined_bams.map { id, files ->
                                           def order = ['N': 0, 'B': 1, 'T': 2]
                                           //tuple(id, *files.sort { order[(it.name =~ /[NBT]/)[0]] })
                                           tuple(id, *files.sort { f ->
                                               //def letter = (f.name =~ /[NBT]/)[0]
                                               def letter = f.name.find(/[NBT]/)
                                               order[letter]
                                            })
                                        }
       //sorted_joined_bams.view()
       idx_filtered_bam = index_filtered_bams(sorted_joined_bams)
       idx_filtered_bam.view()

       // Msisensor pro
       msisensor_pro(idx_filtered_bam,ref_genome_fasta_ch)                                                  
}

workflow {

    // set channels
    target_dir_ch = Channel.fromPath(params.target_dir, type: "dir")
    wes_bedfile_ch = Channel.value(params.wes_bedfile)
    tmpdir_ch = Channel.value(params.tmpdir)
    ref_genome_fasta_ch = Channel.value(params.ref_genome_fasta)
    ref_genome_fai_ch = Channel.value(params.ref_genome_fai)
    ref_genome_dict_ch = Channel.value(params.ref_genome_dict)
    library_type_ch = Channel.value(params.library_type)
    outdir_ch = Channel.value(params.outdir)

    postprocessing_wes_custom_oa(target_dir_ch,wes_bedfile_ch,tmpdir_ch,ref_genome_fasta_ch,ref_genome_fai_ch,ref_genome_dict_ch,library_type_ch,outdir_ch)
}
