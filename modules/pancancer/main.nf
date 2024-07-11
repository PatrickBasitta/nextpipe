nextflow.enable.dsl=2

process ENSEMBL_VEP {
  conda "bioconda::ensembl-vep=111.0"
 
  publishDir "${params.outdir}/${ID}/ENSEMBL_VEP", mode: "copy"
  
  input:
  tuple val(ID), path(vcf)
  val(vep_cache)
  val(fasta)

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
  val(transcript_lst)
  val(variantDBi)

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

workflow PANCANCER_CLC_VEP {
  take:
    args

  main:
   
    def vcfs = args.vcf
    def clc_csvs = args.clc_csv
    def dir_cache = args.vep_cache
    def fasta = args.fasta
    def transcript_list = args.transcriptlist
    def variantDBi = args.variantlist
    def outdir = args.outdir 
    println variantDBi
                 
    vcf_ch = Channel.fromPath(vcfs).map { vcf_file  -> [vcf_file.getSimpleName(), vcf_file]}
     
    dir_cache_ch = Channel.value(dir_cache)
    fasta_ch = Channel.value(fasta)
    
    clc_csv_ch = Channel.fromPath(clc_csvs, checkIfExists: true).map { clc_file -> [clc_file.getSimpleName(), clc_file]}
            
    ENSEMBL_VEP(vcf_ch,dir_cache_ch,fasta_ch)

    format_vep_outputs(ENSEMBL_VEP.out.tab)

    csv_vcf_ch = clc_csv_ch.join(format_vep_outputs.out.formatted_vep_data)
    //variantDBi_ch.view()
    //organize_variant_table(csv_vcf_ch,transcript_lst_ch,args.variantlist)
    organize_variant_table(csv_vcf_ch, transcript_list, args.variantlist)
    final_table = organize_variant_table.out.final_xlsx
  emit:
    final_table
  
 
}

workflow {
  def args = [:]
  for (param in params) { args[param.key] = param.value }
  PANCANCER_CLC_VEP(args)
}

