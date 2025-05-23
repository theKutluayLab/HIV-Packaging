#!/bin/bash

declare -a arr=(
'/path/to/BBduk_folder/' #with slash
) 

index="/path/to/reference_genome_folder"
gtf="/path/to/reference_genome.gtf"
script="/path/to/align_bowtie_ris.sh"

for i in "${arr[@]}"
do
  baseoutput=$(dirname "$i")

  cd "$i" || { echo "cd Failure"; exit 1; }
  find . -maxdepth 1 -name "BC*" -print | while read file
  do
    n=$(sed -e 's#.*BC\(\)#\1#' <<< "$file" | cut -f 1 -d '.')
    fqfile=$i'BC'$n'.fastq'
    basepath=$(basename "$fqfile" .fastq)'-'$(basename $index)
    echo "$basepath"
    output=$baseoutput'/Bowtie_hgenome_mapping/'$(basename "$fqfile" .fastq)'/'

    bash $script $index "$fqfile" "$output" "$basepath" $gtf
  done

done
