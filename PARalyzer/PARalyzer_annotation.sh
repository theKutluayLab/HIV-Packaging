#!/bin/bash

exp_name="Name of the experiment"
input_dir="/path/to/PARalyzer_input_folder"
out_dir="/path/to/PARalyzer_output_folder/${exp_name}"

gtf_dir="/path/to/mart_export1.gtf"
pl_dir="/path/to/folder/PARalyzer-cluster-annotation"

annotated_peaks_folder="$out_dir/AnnotatedPeaks"
annotated_summary_folder="$out_dir/AnnotatedSummary"
mkdir -p "$annotated_peaks_folder"
mkdir -p "$annotated_summary_folder"

cluster_csv="$input_dir/${exp_name}.csv"
output_file="$annotated_peaks_folder/${exp_name}_anno.txt"
output_anno_summary="$annotated_summary_folder/${exp_name}_anno_summary.txt"
output_anno_tRNA="$annotated_summary_folder/${exp_name}_anno_tRNA.txt"

perl "$pl_dir"/Anotbingenereads.pl "$gtf_dir" "$cluster_csv" "$output_file"
perl "$pl_dir"/summarize_annotation_detailed.pl "$output_file" > "$output_anno_summary"
perl "$pl_dir"/summarize_tRNA.pl "$output_file" > "$output_anno_tRNA"