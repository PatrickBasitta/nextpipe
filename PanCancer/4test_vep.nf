nextflow.enable.dsl=2

process VEP_TOOL {

  conda "bioconda::ensembl-vep=108.0"
  
  input:
  path(vcf)
  path(vep_cache)
  //val genome
  //val species
  //val cache_version
  path(fasta)

  output:
  path("*.txt"),  emit: tab
  path("*summary.html"), optional:true, emit: report
  path("*.yml"), optional:true, emit: version

  script:
  """ 
  vep -i ${vcf} -o ${vcf}.txt --tab --everything --species homo_sapiens --assembly GRCh38 \
      --merged --format vcf --force_overwrite --cache_version 108 --cache --dir_cache ${vep_cache} \
      --fasta ${fasta} --offline --dir_plugins ${vep_cache}/plugins/ --plugin SpliceRegion 
  """
}

process FORMAT_VEP_DATA {

  input:
  path(vep_data)

  output:
  path("*.txt"), emit: formatted_vep_data

  script:
  """
  cat ${vep_data} | grep -v '##' > formatted_${vep_data}.txt
  """
}

process PANCANCER_PROCESSING {

  input:
  path(clc_txt)
  path(formatted_vep_data)
  path(transcript_lst)
  path(outfile)

  output:
  path("*.xlsx"), emit: combined_vep_clc
  // stdout
  script:
  """
  python /data2/basitta/bin/PAN_varianten_CLI.py --clc_PAN_file ${clc_txt} --vep_PAN_file ${formatted_vep_data} --transcript_PAN_list ${transcript_lst} --outfile ${outfile}.xlsx
  """
}


workflow VEP {
  take:
    args

  main:
    def vcf = args.vcf
    def dir_cache = args.dir_cache
    def fasta = args.fasta
    def clc_txt = args.clc_txt
    def transcript_lst = args.transcript_lst
    //def outdir = args.outdir
    def outfile = args.outfile

    vcf_ch = Channel.fromPath(vcf)
    dir_cache_ch = Channel.fromPath(dir_cache)
    fasta_ch = Channel.fromPath(fasta)
    clc_txt_ch = Channel.fromPath(clc_txt)
    transcript_lst_ch = Channel.fromPath(transcript_lst)
    //outdir_ch = Channel.fromPath(outdir)
    outfile_ch = Channel.fromPath(outfile)
    
    VEP_TOOL(vcf_ch,dir_cache_ch,fasta_ch)
    FORMAT_VEP_DATA(VEP_TOOL.out.tab)
    PANCANCER_PROCESSING(clc_txt_ch,FORMAT_VEP_DATA.out.formatted_vep_data,transcript_lst_ch,outfile_ch)
    //VEP_TOOL.out.tab.view()
    //vcf_ch.view()
    //dir_cache_ch.view()
    //fasta_ch.view()
    PANCANCER_PROCESSING.out.view() 
}

workflow {
  def args = params
  VEP(args)  
}

