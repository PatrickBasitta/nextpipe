nextflow.enable.dsl=2

process ENSEMBL_VEP {

  tag "${ID}"
  conda "bioconda::ensembl-vep=111.0"
  debug true
 
  publishDir "${params.outdir}/${ID}/ENSEMBL_VEP", mode: "copy"
  
  input:
  tuple val(ID), path(vcf)
  val(vep_cache)
  val(fasta)

  output:
  tuple val(ID), path("*.txt")        ,  emit: tab
  tuple val(ID), path("*summary.html"),  emit: report
        
  script:
  """ 
  mkdir -p ${params.outdir}/${ID}
  vep -i ${vcf} -o ${vcf}.txt --tab --everything --species homo_sapiens --assembly GRCh38 \
      --merged --format vcf --force_overwrite --cache_version 111 --cache --dir_cache ${vep_cache} \
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
  tuple val(ID), path(csv), path(formatted_vep_data)
  val(transcript_lst)
  val(variantDBi)

  output:
  path("*_final_processed.xlsx"), emit: final_xlsx
  path("*_removed_variants.xlsx"), emit: removed_variants
  path("*.log"), emit: log
  
  script:
  """
  python /data/basitta/nextpipe/PanCancer/bin/PAN_varianten_v1.0.py \
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
   
    def vcfs = args.vcfs
    def clc_csvs = args.clc_csvs
    def dir_cache = args.dir_cache
    def fasta = args.fasta
    def transcript_lst = args.transcript_lst
    def variantDBi = args.variantDBi
    def outdir = args.outdir 
                 
    vcf_ch = Channel
                 .fromPath(vcfs)
                 .map { vcf_file  -> [vcf_file.getSimpleName(), vcf_file]}
    vcf_ch.view() 
    dir_cache_ch = Channel.value(dir_cache)
    fasta_ch = Channel.value(fasta)
    
    clc_csv_ch = Channel
                 .fromPath(clc_csvs, checkIfExists: true)
                 .map { clc_file -> [clc_file.getSimpleName(), clc_file]}
    clc_csv_ch.view() 
    transcript_lst_ch = Channel.value(transcript_lst)
    variantDBi_ch = Channel.value(variantDBi)
            
    ENSEMBL_VEP(vcf_ch,dir_cache_ch,fasta_ch)
    FORMAT_VEP_DATA(ENSEMBL_VEP.out.tab)
    csv_vcf_ch = clc_csv_ch.join(FORMAT_VEP_DATA.out.formatted_vep_data)
    PANCANCER_PROCESSING(csv_vcf_ch,transcript_lst_ch,variantDBi_ch)
 
}

workflow {
  def args = params
  PANCANCER_CLC_VEP(args)
}

