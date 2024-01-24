nextflow.enable.dsl=2

process ENSEMBL_VEP {

  tag "${ID}"
  conda "bioconda::ensembl-vep=108.0"
  debug true
 
  publishDir "${params.outdir}/${ID}/ENSEMBL_VEP", mode: "copy"
  
  input:
  tuple val(ID), path(vcf)
  path(vep_cache)
  path(fasta)

  output:
  tuple val(ID), path("*.txt")        ,  emit: tab
  tuple val(ID), path("*summary.html"),  emit: report
      
  script:
  """ 
  mkdir -p ${params.outdir}/${ID}
  vep -i ${vcf} -o ${vcf}.txt --tab --everything --species homo_sapiens --assembly GRCh38 \
      --merged --format vcf --force_overwrite --cache_version 108 --cache --dir_cache ${vep_cache} \
      --fasta ${fasta} --offline --dir_plugins ${vep_cache}/plugins/ --plugin SpliceRegion 
  """
}

process FORMAT_VEP_DATA {
  
  tag "${ID}"
  debug true

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

process PANCANCER_PROCESSING {
  
  tag "${ID}"
  debug true
 
  publishDir "${params.outdir}/${ID}/FINAL_OUTPUT", mode: "copy"
  
  input:
  tuple val (ID), path(csv), path(formatted_vep_data)
  path(transcript_lst)

  output:
  path("*.xlsx"), emit: final_xlsx
  
  script:
  """
  python /data2/basitta/bin/PAN_varianten_CLI.py --clc_PAN_file ${csv} --vep_PAN_file ${formatted_vep_data} --transcript_PAN_list ${transcript_lst} --outfile ${ID}_final_processed.xlsx
  """
}


workflow PANCANCER_CLC_VEP {
 
  take:
    args

  main:
    def vcf = args.vcf
    def dir_cache = args.dir_cache
    def fasta = args.fasta
    def clc_txt = args.clc_txt
    def transcript_lst = args.transcript_lst
    def outdir = args.outdir
   
    vcf_ch = Channel
                 .fromPath(vcf, checkIfExists: true)
                 .map { vcf_file -> [vcf_file.getSimpleName(), vcf_file]}
    dir_cache_ch = Channel.fromPath(dir_cache)
    fasta_ch = Channel.fromPath(fasta)
    clc_csv_ch = Channel
                 .fromPath(clc_txt, checkIfExists: true)
                 .map { clc_file -> [clc_file.getSimpleName(), clc_file]}
    transcript_lst_ch = Channel.fromPath(transcript_lst)
    
    vcf_ch.view()
    clc_csv_ch.view()   

    ENSEMBL_VEP(vcf_ch,dir_cache_ch,fasta_ch)
    FORMAT_VEP_DATA(ENSEMBL_VEP.out.tab)
    csv_vcf_ch = clc_csv_ch.join(FORMAT_VEP_DATA.out.formatted_vep_data)
    PANCANCER_PROCESSING(csv_vcf_ch,transcript_lst_ch)
    PANCANCER_PROCESSING.out.view() 
}

workflow {
  def args = params
  PANCANCER_CLC_VEP(args)  
}

