#!/bin/bash
prefix="Variants-"
target_dir_csv="/PAT-Sequenzer/PANCANCER/INPUT/csv_input"
target_dir_vcf="/PAT-Sequenzer/PANCANCER/INPUT/vcf_input"

mv -v /PAT-Sequenzer/ImportExport/CLC_PAN_EXPORT/*.vcf $target_dir_vcf
mv -v /PAT-Sequenzer/ImportExport/CLC_PAN_EXPORT/*.csv $target_dir_csv

for filevcf in /PAT-Sequenzer/PANCANCER/INPUT/vcf_input/*.vcf
do
    SUBSTRING=$(echo $filevcf | cut -d' ' -f4 | cut -d'_' -f1)
    ID_catched=${SUBSTRING#"$prefix"}
    mv -v "$filevcf" ${target_dir_vcf}/"$ID_catched.vcf"
done

for filecsv in /PAT-Sequenzer/PANCANCER/INPUT/csv_input/*.csv
do
    SUBSTRING=$(echo $filecsv | cut -d' ' -f4 | cut -d'_' -f1)
    ID_catched=${SUBSTRING#"$prefix"}
    mv -v "$filecsv" ${target_dir_csv}/"$ID_catched.csv"
done
