#!/bin/bash

infile=$1 # the fq input file
barcode=$2 
fd=$3
adapter=$4
#cd /opt/conda/bin || { echo "cd Failure"; exit 1; }


tp=$fd'temp/'
logdir=$fd'logs/'
miss=$fd'miss_filter/'

mkdir -p "$tp"
mkdir -p "$logdir"
mkdir -p "$miss"

#The code below reads barcodes from a barcode file into an array

declare -a indnum=()
declare -a bcnum=()
declare -a codes=()
readarray -t barcodes < "$barcode"

for i in "${!barcodes[@]}"
do
	line=(${barcodes[i]//;/ }) 
	indnum[i]=$i
	bcnum[i]=${line[0]}
	codes[i]=${line[1]}
done

# collapsing identical reads
dedupe.sh in="$infile" out="$tp"'collap.fq' ac=f 2> "$logdir"'dedupe_log.txt'

# Trimming adapters
bbduk.sh -Xmx51g in="$tp"'collap.fq' out="$tp"'atrim.fq' outm="$tp"'atrim_miss.fq' literal="$adapter" ktrim=r k=8 rcomp=f mink=7 ml=18 maxlength=52 2> "$logdir"'ATlog.txt'


for i in "${indnum[@]}"
do
	bbduk.sh -Xmx51g in="$tp"'collap.fq' outm="$tp""${bcnum[$i]}"'t.fq' literal='NN'"${codes[$i]}" k=8 rcomp=f copyundefined mm=f 2> "$logdir"'filter'"${bcnum[$i]}".txt
	bbduk.sh -Xmx51g in="$tp""${bcnum[$i]}"'t.fq' out="$fd""${bcnum[$i]}".fq outm="$miss""${bcnum[i]}"'_miss.fq' literal='NN'"${codes[$i]}" copyundefined ktrim=r k=8 rcomp=f mm=f ml=10 maxlength=44 2> "$logdir"'trim'"${bcnum[$i]}".txt
done


cd "$logdir" || { echo "cd Failure"; exit 1; }
cat * > merged-file.txt
