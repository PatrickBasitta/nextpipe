#!/bin/bash

nextflow run PANCANCER.nf \
         -c resources.config \
         --vcfs  '/PAT-Sequenzer/PANCANCER_INPUT_NXF/vcf_input/*.vcf' \
         --clc_csvs '/PAT-Sequenzer/PANCANCER_INPUT_NXF/csv_input/*.csv' \
         --dir_cache /data2/basitta/vep_cache/ \
         --fasta /data2/basitta/vep_fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
         --transcript_lst /data2/basitta/pan_data/finale_aktuelle_pancancer_transcript_lst.xlsx \
         --variantDBi /PAT-Sequenzer/PanCancer_test/Variantenliste22_12_15.xlsx \
         --outdir /PAT-Sequenzer/PANCANCER_OUTPUT_NXF/ \
         -with-conda \
         -resume

