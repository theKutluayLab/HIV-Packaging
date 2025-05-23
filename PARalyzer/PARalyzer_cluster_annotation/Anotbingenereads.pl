#! /usr/bin/perl
##############################################################################################
# bingenereads v1.0 Anjan Purkayastha                                                    #
# A program to bin reads in a .rlf file for each gene. Gene information entered             # 
#according to the refseq format.                                                            #
#Usage: perl binGeneReads.pl read_file.rlf gene_annotation_file.txt > output_gene_reads file#
#############################################################################################

use strict;
use warnings;

open (R, "$ARGV[0]") || die "cannot open file to read: $!";

##############################################################################################
#Create a Hash of hashes by parsing the information                                          # 
#from the Annotation files.                                                                  #
#GnID GnName TrID PrID Chr Strand 5UTRst 5UTRend 3UTRst 3UTRend ExStart ExEnd ExRank Biotype #
#	0	1		2	3	4	5		6		7		8		9	10		11		12  	13   #
##############################################################################################
my %HOH;
my $header_skip= <R>;
my $i =0;
my $b = 0;
my $TrID = '';
my $baseIn = '';
my $prev;
my $prevEx;
while(<R>){
	$i++;
	chomp;
	my @line= split /\t/;
	my $chr= "chr" . $line[4]; #add "chr" to chromsome name, chromosomes named as chrXX in the Refseq file. However in the .RLF format only chromosome number is stored.
	#my $chr= $line[0]; 
	my $base5UTR = '';
	my $base3UTR = '';
	my $baseEx;
	my $exst = $line[10];
	my $exend = $line[11];
	
	if ($line[6]) {
		if ($line[5] eq "-1"){	
			$base5UTR = "-" . $line[6] . '-' . $line[7];
			$exend = $line[6] - 1;
		}else{
			$base5UTR = $line[6] . '-' . $line[7];
			$exst = $line[7] + 1;
		}
	}
	
	if ($line[8]) {
		if ($line[5] eq "-1"){	
			$base3UTR = "-" . $line[8] . '-' . $line[9];
			$exst = $line[9] + 1;
		}else{
			$base3UTR = $line[8] . '-' . $line[9];
			$exend = $line[8] - 1;
		}
	}
	
	if ($line[5] eq "-1"){	
		$baseEx = "-" . $exst . '-' . $exend;
	}else{
		$baseEx = $exst . '-' . $exend;
	}
	my $anot = $line[1] . '|' . $line[13] . '|' . $line[2] . '|';
	
	my $currTrID = $line[2];
	
	if($base5UTR){
		my $anotUTR = $anot . '5UTR';
		$HOH{$chr}{$base5UTR} = $anotUTR;
	}elsif($base3UTR){
		my $anotUTR = $anot . '3UTR';
		$HOH{$chr}{$base3UTR} = $anotUTR;
	}

	if($TrID ne $currTrID){
		$TrID = $currTrID;
		if ($line[5] eq "-1"){
			$prev = $line[10];
		}else{
			$prev = $line[11];
		}
		$prevEx = $line[12];
	}else{
		$TrID = $currTrID;
		if ($line[5] eq "-1"){
			$baseIn = "-" . $line[11] . '-' . $prev;
			$prev = $line[10];
		}else{
			$baseIn = $prev . '-' . $line[10];
			$prev = $line[11];
		}
		$HOH{$chr}{$baseIn} = $anot . 'Intron';
		$prevEx = $line[12];
	}

	
	$anot .= 'Exon';
	
	
	#$total+=$count;
	$HOH{$chr}{$baseEx}= $anot;
	
	if($i > 100){
		$i = 1;
		$b++;
		my $printnum = ($b * 100);
		print "processed $printnum\n";
	}
}


#my $norm_factor= $total/(10**6); #normalization factor to divide number of reads by the number of million mapped reads.

#############################################################
#Parse the custers list.                                    #
#Count number of reads found within the coordinates of each #
#gene. Print Gene name and normalized read count.           #
#############################################################
open (G, "$ARGV[1]") || die "cannot open refseq file to read: $!";
open (OUT, ">$ARGV[2]") || die "cannot make results!!";
print OUT "Chromosome\tStrand\tClusterStart\tClusterEnd\tClusterID\tClusterSequence\tReadCount\tModeLocation\tModeScore\tConversionLocationCount\tConversionEventCount\tNonConversionEventCount\tGene\tBiotype\tTrID\tUTR/Exon\n";
$i =0;
$b =0;
my $header_skip2= <G>;
while (<G>){
	chomp;
	$i++;
	my @line= split/,/;
	my $linet = join("\t",@line);
	my $target_chr= $line[0];
	#my $gene= $line[0];
	#my $id= $line[1];
	my $strand= $line[1];
	my $left= $line[2];
	my $right=$line[3];
	my $halfseqlen = int(length($line[5]) / 2);
	#my $read_count=0;
	my %priority = ();
	if ($strand eq "+"){	
		for my $ref (keys %{$HOH{$target_chr}}){
			if ($ref !~ /^-/){
				my @coo = split(/-/,$ref);
				if ((($coo[0] - $halfseqlen) <= $left) && (($coo[1] + $halfseqlen) >= $right)){
					my @anotat = split(/\|/,$HOH{$target_chr}{$ref});
					$priority{$anotat[3]} = "$linet\t$anotat[0]\t$anotat[1]\t$anotat[2]\t$anotat[3]\n";
					#print ("$ref\t$HOH{$target_chr}{$ref}\n");
					#$read_count+= $HOH{$target_chr}{$ref};
				}
			}
		}
	}

	if ($strand eq "-"){
		for my $ref (keys %{$HOH{$target_chr}}){
			if($ref =~ /^-/){
				#my $original_ref= $ref;
				#$ref=~ s/-//;
				my @coo = split(/-/,$ref);
				if ((($coo[1] - $halfseqlen) <= $left) && (($coo[2] + $halfseqlen) >= $right)){
					my @anotat = split(/\|/,$HOH{$target_chr}{$ref});
					$priority{$anotat[3]} = "$linet\t$anotat[0]\t$anotat[1]\t$anotat[2]\t$anotat[3]\n";
					#print ("$original_ref\t$HOH{$target_chr}{$original_ref}\n");
					#$read_count+= $HOH{$target_chr}{$original_ref};
				}
			}
		}
	}

#if ($read_count ==0){ # set a floor of total read count at 1.This will avoid difficulties while working in log space.
	#$read_count=1;
#}

#my $norm_read= $read_count/$norm_factor;
#my $norm_read_final= sprintf("%10.2f", $norm_read);

#print ("$gene\t$id\t$strand\t$target_chr\t$read_count\n");
	if ($priority{'Exon'}){
		print OUT "$priority{'Exon'}";
	}elsif ($priority{'5UTR'}){
		print OUT "$priority{'5UTR'}";
	}elsif ($priority{'3UTR'}){
		print OUT "$priority{'3UTR'}";
	}elsif ($priority{'Intron'}){
		print OUT "$priority{'Intron'}";
	}else{
		print OUT "$linet\tNA\tNA\tNA\tIntergenic\n";
	}
	
	if($i > 100){
		$i = 1;
		$b++;
		my $printnum = ($b * 100);
		print "cluster $printnum\n";
	}


}


