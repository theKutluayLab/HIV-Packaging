#!/bin/bash

#BBDuk file path
declare -a arr=('/path/to/.fastq')

#modify the directory below to change barcodes
barcode_file="/storage1/fs1/kutluay/Active/home/Inputs/barcodes_maritza.txt"
script="/storage1/fs1/kutluay/Active/home/Scripts_10_2019/project_scripts/KV-Exp250-252-CLIP-seq-2-C/mapping/BBduk/bbduk_run_rp_ris_10nt_collapse.sh"

for i in "${arr[@]}"
do
  base=$(basename "${i%.*}")
  output='/path/to/output_folder'

  mkdir -p "$output"
  adapter='TGGAATTC'
  bash $script $i $barcode_file $output $adapter
done