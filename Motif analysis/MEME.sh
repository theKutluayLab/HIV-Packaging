#!/bin/bash

# input Parameters
input_cluster="/path/to/-cluster-Virion.txt"
output_bed="/path/to/output.bed"
gtf_file="/path/to/reference_genome.gtf"
output_mrna_gtf="/path/to/mRNA_gtf_output.gtf"
bed_dir="/path/to/bed_files"
ref_genome="/path/to/reference_genome.fa"
fa_output_dir="/path/to/output_fa_folder"
meme_output_base="/path/to/meme_output_folder"


# convert PARalyzer output to BED
paralyzer_to_bed() {
  local input_file="$1"
  local output_file="$2"
  awk 'BEGIN {OFS="\t"} NR > 1 && $14 == "protein_coding" {print $1, $3, $4, $6, $2}' "$input_file" > "$output_file"
}

# extract protein-coding mRNA regions from GTF
extract_mrna_from_gtf() {
  local gtf_file="$1"
  local output_file="$2"
  awk '$3 == "transcript" && /gene_biotype "protein_coding"/' "$gtf_file" > "$output_file"
}

# convert BED to deduplicated FASTA
bed_to_fasta_dedup() {
  local bed_dir="$1"
  local genome_fasta="$2"
  local output_dir="$3"
  mkdir -p "$output_dir"

  for bed_file in "$bed_dir"/*.bed; do
    local base_name
    base_name=$(basename "$bed_file" .bed)
    local temp_fa="$output_dir/$base_name.fa"
    local final_fa="$output_dir/${base_name}_dedup.fa"
    bedtools getfasta -fi "$genome_fasta" -bed "$bed_file" -fo "$temp_fa"
    seqkit rmdup -s "$temp_fa" > "$final_fa"
    rm "$temp_fa"
  done
}

# run MEME
run_meme() {
  local fa_dir="$1"
  local meme_base_output="$2"
  mkdir -p "$meme_base_output"

  for fasta_file in "$fa_dir"/*.fa; do
    local base_name
    base_name=$(basename "$fasta_file" .fa)
    local meme_output_dir="$meme_base_output/$base_name"
    mkdir -p "$meme_output_dir"
    meme "$fasta_file" -oc "$meme_output_dir" -nmotifs 10 -minw 5 -maxw 20 -dna -mod anr
  done
}


#run file
paralyzer_to_bed "$input_cluster" "$output_bed"
extract_mrna_from_gtf "$gtf_file" "$output_mrna_gtf"
bed_to_fasta_dedup "$bed_dir" "$ref_genome" "$fa_output_dir"
run_meme "$fa_output_dir" "$meme_output_base"
