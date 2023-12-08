nextflow.enable.dsl=2
/*
 * Input files
 */
params.vcf_file = "/PAT-Sequenzer/WES-Pilot/WES_Pilot_final/ready_for_submit_draft/WES_Pilot_vcf_export_15_11_2023/vcfs/*.vcf"
//params.qc_csv = ""
params.chainfile = ""
params.refgenome = ""
params.outdir = ""
params.cpus = 12
params.memory = 60

/*
 * Channels
 */

//qc_csv_ch = Channel.fromPath(params.qc_csv) | splitCsv(header: true) | map { row -> tuple(row.sampleid,row.bamfile) }
samples_ch = Channel.fromPath(params.vcf_file) | map { tuple(it) }
chainfile_ch = Channel.fromPath(params.chainfile)
refgenome_ch = Channel.fromPath(params.refgenome)

process WES_CrossMap {
    
    tag "${vcf_file.getSimpleName()}"
    cpus "${params.cpus}"
    memory "${params.memory}"     
    debug true

    publishDir "${params.outdir}", mode: "copy"      
    
    input:
        path(vcf_file)
        path(chainfile)
        path(refgenome)
    output:
        path("${vcf_file.getSimpleName()}_remapped.vcf")
    script:
    """
    CrossMap.py vcf ${chainfile} ${vcf_file} ${refgenome} ${vcf_file.getSimpleName()}_remapped.vcf 
    """
} 

//workflow WES_VCF_REMAPPED {
//    take:
//        samples_ch
//        chainfile_ch
//        refgenome_ch  
//    main:
//        WES_CrossMap(samples_ch,chainfile_ch,refgenome_ch)
//    emit:
//        remapped = WES_CrossMap.out
//}
    
//workflow {
    //samples_ch.view()
//    WES_VCF_REMAPPED(samples_ch,chainfile_ch,refgenome_ch)
//}

workflow {
    WES_CrossMap(samples_ch,chainfile_ch,refgenome_ch)
}
