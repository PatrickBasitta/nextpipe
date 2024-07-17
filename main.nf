nextflow.enable.dsl=2

process ensembl_vep {
  cpus 10
  memory 35.GB

  conda "bioconda::ensembl-vep=111.0"
 
  publishDir "${params.outdir}/${ID}/ENSEMBL_VEP", mode: "copy"

  input:
  tuple val(ID), path(vcf)
  path(vep_cache)
  path(fasta)

  output:
  //tuple val(ID), path("*.txt")        ,  emit: tab
  tuple val(ID), path("${vcf}.txt")        ,  emit: tab
  tuple val(ID), path("*summary.html"),  emit: report
        
  script:
  """ 
  mkdir -p ${params.outdir}/${ID}
  vep -i ${vcf} -o ${vcf}.txt --tab --everything --species homo_sapiens --assembly GRCh38 \
      --merged --format vcf --force_overwrite --cache_version 111 --cache --dir_cache ${vep_cache} \
      --fasta ${fasta} --offline --dir_plugins ${vep_cache}/plugins/ --plugin SpliceRegion 
  """
}

process format_vep_outputs {
  publishDir "${params.outdir}/${ID}/FORMAT_VEP_DATA", mode: "copy"  

  input:
  tuple val(ID), path(vep_data)

  output:
  tuple val(ID), path("*formatted_*"), emit: formatted_vep_data
  tuple val(ID), path("vcf_header_*")

  script:
  """
  cat ${vep_data} | grep -v '##' > formatted_${vep_data}
  cat ${vep_data} | grep '##' > vcf_header_${vep_data}
  """
}

// format, filter, merge and unify data to make it readable
// and interpretable when looking at the table
process organize_variant_table {
  conda "pandas=2.2.2 openpyxl=3.1.4"
 
  publishDir "${params.outdir}/${ID}/FINAL_OUTPUT", mode: "copy"
  
  input:
  tuple val(ID), path(csv), path(formatted_vep_data)
  path(transcript_lst)
  path(variantDBi)

  output:
  path("${ID}_final_processed.xlsx"), emit: final_xlsx
  path("${ID}_removed_variants.xlsx"), emit: removed_variants
  path("log_${ID}.log"), emit: log
  
  script:
  """
  organize_variant_table.py \
      --clc ${csv} \
      --vep ${formatted_vep_data} \
      -t ${transcript_lst} \
      -D ${variantDBi} \
      -o ${ID}_final_processed.xlsx \
      -rv ${ID}_removed_variants.xlsx > log_${ID}.log
  """
}

workflow pancancer_analyse {
  take:
    vcfs
    clc_csvs
    dir_cache
    fasta
    transcript_list
    variantlist
    outdir
  main:
    vcf_ch = vcfs.map { vcf_file  -> [vcf_file.getSimpleName(), vcf_file]}
     
    dir_cache_ch = dir_cache
    fasta_ch = fasta
    
    clc_csv_ch = clc_csvs.map { clc_file -> [clc_file.getSimpleName(), clc_file]}
            
    tab = ensembl_vep(vcf_ch,dir_cache_ch,fasta_ch).tab

    format_vep_outputs(tab)

    csv_vcf_ch = clc_csv_ch.join(format_vep_outputs.out.formatted_vep_data)
    organize_variant_table(csv_vcf_ch, transcript_list, variantlist)
    final_table = organize_variant_table.out.final_xlsx
  emit:
    final_table
  
 
}

workflow {
  def args = [:]
  for (param in params) { args[param.key] = param.value }
  pancancer_analyse(
  	args.vcf, args.clc_csv, args.vep_cache, args.fasta, 
  	args.transcriptlist, args.variantlist, args.outdir
	)
}

