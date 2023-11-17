# Process/Algorithm processing exome data for WES-Pilot

# Step1: Fastq to vcfs
#        Run CLC Exome Workflow (custom Qiaqen workflow) with fastq WES data
#        Output: VariantTrack (export as vcf; genome build: hg38) 

# Step2: Upload vcf to QCI, submit, interpret and download corresponding vcf,tsv,csv file

# Step3: Upload QCI vcf file to CLC Workbench

# Step4: Download CLC tracks - QCI and CLC Variant track - as txt.files in \\NAS: WES-Pilot/WES_Pilot_input_all (in Windows)

# Step5: Run WES_Pilot_hg38_CLI.py and interpret final xlsx file

# Step6: NXF script to report in WES format (here liftover hg38 to hg19)

