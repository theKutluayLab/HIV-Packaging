#!/bin/bash

index=$1
fqfile=$2
output=$3
basepath=$4
gtf=$5

mkdir -p "$output"
cd "$output" || { echo "cd Failure"; exit 1; }

#bowtie "$index" -q "$fqfile" -v 1 -m 10 --best --strata -S "$output$basepath"'.sam'  # report multimapping reads
bowtie "$index" -q "$fqfile" -v 1 -m 10 --best --strata -S --norc "$output$basepath"'.sam' # only map to sesne strand


samtools view -b -h -o "$output$basepath"'.bam' "$output$basepath"'.sam'
samtools sort -o "$output$basepath"'.sorted.bam' "$output$basepath"'.bam'
samtools index "$output$basepath"'.sorted.bam'
samtools view -h -F 4 "$output$basepath"'.sorted.bam' > "$output$basepath"'.sorted.sam' #newly added for annotating peaks in PARalyzer
samtools view -b -f 4 "$output$basepath"'.sorted.bam' > "$output"'viral_unmapped.bam'
samtools fastq "$output"'viral_unmapped.bam' > "$output"'viral_unmapped.fq'

#generate .bed files 
bed_file="$output$basepath.sorted.bed"
bedtools bamtobed -i "${output}${basepath}.sorted.bam" > "$bed_file"

featureCounts -t exon -g gene_id -a "$gtf" -o "$output$basepath"'_fcCounts.count' "$output$basepath"'Aligned.sortedByCoord.out.bam'
python /path/to/reduce_fc.py "$output$basepath"'_fcCounts.count'