#!/bin/bash
archive_dir_csv="/PAT-Sequenzer/PANCANCER/ARCHIVE_AND_LOGS/csv_archiv"
archive_dir_vcf="/PAT-Sequenzer/PANCANCER/ARCHIVE_AND_LOGS/vcf_archiv"

mv -v /PAT-Sequenzer/PANCANCER/INPUT/vcf_input/*.vcf $archive_dir_vcf
mv -v /PAT-Sequenzer/PANCANCER/INPUT/csv_input/*.csv $archive_dir_csv
