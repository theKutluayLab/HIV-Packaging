#!/usr/bin/perl

use warnings;
use strict;

my $file = $ARGV[0] or die "Need to get a tab-delimited file on the command line\n";
open my $fh, '<', $file;
my $intergenic=0;
my $rintergenic=0;
my $genic=0;
my $rgenic=0;
my $exon=0;
my $rexon=0;
my $intron=0;
my $rintron=0;
my $tRNA=0;
my $rtRNA=0;
my $miscRNA=0;
my $rmiscRNA=0;
my $rRNA=0;
my $rrRNA=0;
my $CDS=0;
my $rCDS=0;
my $FiveUTR=0;
my $r5UTR=0;
my $ThreeUTR=0;
my $r3UTR=0;
my $SevenSL=0;
my $r7SL=0;
my $mRNA=0;
my $rmRNA=0;
my $RN7SL1=0; 
my $rRN7SL1=0;
my $RN7SL2=0;
my $rRN7SL2=0;
my $RN7SL3=0;
my $rRN7SL3=0;
my $miRNA=0;
my $rmiRNA=0;
my $lincRNA=0;
my $rlincRNA=0;

while (<$fh>) {

	chomp($_);
	my @line = split(/\t/,$_);
	
if ($line[15] eq 'Intergenic'){
        $intergenic++;
	$rintergenic += $line[6];
        }
if ($line[15] eq '5UTR'){
	$FiveUTR++;
	$r5UTR += $line[6];
	}
if ($line[15] eq '3UTR'){
        $ThreeUTR++;
	$r3UTR += $line[6];
        }
if ($line[15] eq 'Exon'){
        $exon++;
        $rexon += $line[6];

	if ($line[13] =~ /tRNA/){
        	$tRNA++;
        	$rtRNA += $line[6];
        }
	if ($line[13] =~ /misc_RNA/){
        	$miscRNA++;
        	$rmiscRNA += $line[6];
        }
	if ($line[12] =~ /RN7SL/){
        	$SevenSL++;
        	$r7SL += $line[6];
		if ($line[13] =~ /antisense/){
                	$miscRNA++;
                	$rmiscRNA += $line[6];
                }
	}
	if ($line[12] eq 'RN7SL1'){
        	$RN7SL1++;
        	$rRN7SL1 += $line[6];
        }
	if ($line[12] eq 'RN7SL2'){
        	$RN7SL2++;
        	$rRN7SL2 += $line[6];
        }
	if ($line[12] eq 'RN7SL3'){
        	$RN7SL3++;
       		$rRN7SL3 += $line[6];
        }
	if ($line[13] =~ /miRNA/){
        	$miRNA++;
        	$rmiRNA += $line[6];
        }
	if ($line[13] =~ /lincRNA/){
        	$lincRNA++;
        	$rlincRNA += $line[6];
        }
	if ($line[13] =~ /rRNA/){
        	$rRNA++;
        	$rrRNA += $line[6];
        }
        if ($line[13] =~ /protein_coding/){
        	$CDS++;
        	$rCDS += $line[6];
	}                                                                                                                                                                        
}

if ($line[15] eq 'Intron'){
        $intron++;
        $rintron += $line[6];
        }
	
}

$genic = $FiveUTR+$ThreeUTR+$intron+$exon;
$rgenic = $r5UTR+$r3UTR+$rintron+$rexon;
print "CLUSTERS:\n";
print "intergenic\t$intergenic\n";
print "genic: 5UTR+3UTR+intron+exon\t$genic\n";
print "5UTR\t$FiveUTR\n";
print "3UTR\t$ThreeUTR\n";
print "Exon\t$exon\n";
print "CDS\t$CDS\n";
print "mRNA: 5UTR+3UTR+CDS\t", $FiveUTR + $ThreeUTR + $CDS, "\n";
print "tRNA\t$tRNA\n";
print "7SL RNA\t$SevenSL\n";
print "rRNA\t$rRNA\n";
print "miRNA\t$miRNA\n";
print "lincRNA\t$lincRNA\n";
print "other_RNA: exon-(tRNA+7SL RNA+miRNA+rRNA+CNS+lincRNA)\t", $exon - ($tRNA+$SevenSL+$miRNA+$rRNA+$CDS+$lincRNA),"\n";
print "Intron\t$intron\n";
print "misc_RNA\t$miscRNA\n";
print "7SL1 RNA\t$RN7SL1\n";
print "7SL2 RNA\t$RN7SL2\n";
print "7SL3 RNA\t$RN7SL3\n";
print "READS WITHIN CLUSTERS:\n";
print "intergenic\t$rintergenic\n";
print "genic: 5UTR+3UTR+intron+exon\t$rgenic\n";
print "5UTR\t$r5UTR\n";
print "3UTR\t$r3UTR\n";
print "Exon\t$rexon\n"; 
print "CDS\t$rCDS\n";
print "mRNA: 5UTR+3UTR+CDS\t", $r5UTR + $r3UTR + $rCDS, "\n";
print "tRNA\t$rtRNA\n";
print "7SLRNA\t$r7SL\n";
print "rRNA\t$rrRNA\n";
print "miRNA\t$rmiRNA\n";
print "lincRNA\t$rlincRNA\n";
print "other_RNA: exon-(tRNA+7SL RNA+miRNA+rRNA+CNS+lincRNA)\t", $rexon - ($rtRNA+$r7SL+$rmiRNA+$rrRNA+$rCDS+$rlincRNA),"\n";
print "Intron\t$rintron\n";
print "misc_RNA\t$rmiscRNA\n";
print "7SL1_RNA\t$rRN7SL1\n";
print "7SL2_RNA\t$rRN7SL2\n";
print "7SL3_RNA\t$rRN7SL3\n";

close $fh;
