#!/bin/bash

index=$1
fqfile=$2
output=$3
basepath=$4
gtf=$5

if [ ! -d "$output" ]
then
mkdir -p "$output"
cd "$output" || { echo "cd Failure"; exit 1; }

bowtie "$index" -q "$fqfile" -v 1 -m 1 -S "$output""$basepath"'.sam' # report uniquely mapped alignments, allow one mismatch
# bowtie "$index" -q "$fqfile" -v 1 -m 10 -S "$output""$basepath"'.sam' # report multi-mapping alignments
# bowtie "$index" -q "$fqfile" -v 2 -m 1 -S "$output""$basepath"'.sam' # report uniquely mapped alignments, allow two mismatches

samtools view -bS "$output""$basepath"'.sam' > "$output""$basepath"'.bam'
samtools sort "$output""$basepath"'.bam' -o "$output""$basepath"'.sorted.bam'
samtools index "$output""$basepath"'.sorted.bam' "$output""$basepath"'.sorted.bai'
echo "Running: samtools mpileup -BQ0 -d 500000 -f $index.fasta $output$basepath.sorted.bam > $output$basepath_mpileup"
samtools mpileup -BQ0 -d 500000 -f "$index"'.fasta' "$output""$basepath"'.sorted.bam' > "$output""$basepath"'_mpileup'

samtools view -b -f 4 "$output""$basepath"'.sorted.bam' > "$output"'viral_unmapped.bam'
samtools fastq "$output"'viral_unmapped.bam' > "$output"'viral_unmapped.fastq'

bed_file="$output""$basepath"'.sorted.bed'
bedtools bamtobed -i "${output}${basepath}.sorted.bam" > "$bed_file"

python3 /storage1/fs1/kutluay/Active/home/Scripts_10_2019/mpileup_to_counts_new.py "$output""$basepath"'_mpileup'
python3 /storage1/fs1/kutluay/Active/home/Scripts_10_2019/pileup_to_counts.py "$output""$basepath"'_mpileup'

for i in $(find . -name '*e*counts*.txt' -exec readlink -f {} \;)
do
  base_name=$(basename "$i" .txt)
  perl /storage1/fs1/kutluay/Active/home/Scripts_10_2019/repair_counts.pl "$i" > "${base_name}_repaired.txt"
done

fi