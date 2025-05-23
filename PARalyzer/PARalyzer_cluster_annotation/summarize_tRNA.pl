#!/usr/bin/perl

use warnings;
use strict;

my $file = $ARGV[0] or die "Need to get a tab-delimited file on the command line\n";
open my $fh, '<', $file;

my $tRNA=0;
my $rtRNA=0;
my $rAlaAGC=0;
my $rAlaCGC=0;
my $rAlaTGC=0;
my $rArgACG=0;
my $rArgCCG=0;
my $rArgCCT=0;
my $rArgTCG=0;
my $rArgTCT=0;
my $rAsnATT=0;
my $rAsnGTT=0;
my $rAspGTC=0;
my $rCysGCA=0;
my $rGlnCTG=0;
my $rGlnTTG=0;
my $rGluCTC=0;
my $rGluTTC=0;
my $rGlyCCC=0;
my $rGlyGCC=0;
my $rGlyTCC=0;
my $rHisGTG=0;
my $rIleAAT=0;
my $rIleGAT=0;
my $rIleTAT=0;
my $rLeuAAG=0;
my $rLeuCAA=0;
my $rLeuCAG=0;
my $rLeuTAA=0;
my $rLeuTAG=0;
my $rLysCTT=0;
my $rLysTTT=0;
my $rMetCAT=0;
my $rPheGAA=0;
my $rProAGG=0;
my $rProCGG=0;
my $rProTGG=0;
my $rPseudo=0;
my $rSeC=0;
my $rSerAGA=0;
my $rSerCGA=0;
my $rSerTGA=0;
my $rSerGCT=0;
my $rSup=0;
my $rThrAGT=0;
my $rThrCGT=0;
my $rThrTGT=0;
my $rTrpCCA=0;
my $rTyrGTA=0;
my $rValAAC=0;
my $rValCAC=0;
my $rValTAC=0;
my $rundet=0;

while (<$fh>) {

	chomp($_);
	my @line = split(/\t/,$_);
	
if ($line[15] eq 'Exon' && $line[13] =~ /tRNA/){
        	 $tRNA++;
                $rtRNA += $line[6];

	if ($line[12] =~ /Ala\(AGC\)/){
   
        	$rAlaAGC += $line[6];
        }
	if ($line[12] =~ /Ala\(CGC\)/){

                $rAlaCGC += $line[6];
	}
	if ($line[12] =~ /Ala\(TGC\)/){

                $rAlaTGC += $line[6];
	}
	if ($line[12] =~ /Arg\(ACG\)/){

                $rArgACG += $line[6];
	}
	if ($line[12] =~ /Arg\(CCG\)/){

                $rArgCCG += $line[6];
        }
	if ($line[12] =~ /Arg\(CCT\)/){

                $rArgCCT += $line[6];
        }
	if ($line[12] =~ /Arg\(TCG\)/){

                $rArgTCG += $line[6];
        }
	if ($line[12] =~ /Arg\(TCT\)/){

                $rArgTCT += $line[6];
        }
	if ($line[12] =~ /Asn\(ATT\)/){

                $rAsnATT += $line[6];
        }
	if ($line[12] =~ /Asn\(GTT\)/){

                $rAsnGTT += $line[6];
        }
	if ($line[12] =~ /Asp\(GTC\)/){

                $rAspGTC += $line[6];
        }
	if ($line[12] =~ /Cys\(GCA\)/){

                $rCysGCA += $line[6];
        }
	if ($line[12] =~ /Gln\(CTG\)/){

                $rGlnCTG += $line[6];
        }
	if ($line[12] =~ /Gln\(TTG\)/){

                $rGlnTTG += $line[6];
        }
	if ($line[12] =~ /Glu\(CTC\)/){

                $rGluCTC += $line[6];
        }
	if ($line[12] =~ /Glu\(TTC\)/){

                $rGluTTC += $line[6];
        }
	if ($line[12] =~ /Gly\(CCC\)/){

                $rGlyCCC += $line[6];
        }
	if ($line[12] =~ /Gly\(GCC\)/){

                $rGlyGCC += $line[6];
        }
	if ($line[12] =~ /Gly\(TCC\)/){

                $rGlyTCC += $line[6];
        }
	if ($line[12] =~ /His\(GTG\)/){

                $rHisGTG += $line[6];
        }
	if ($line[12] =~ /Ile\(AAT\)/){

                $rIleAAT += $line[6];
        }
	if ($line[12] =~ /Ile\(GAT\)/){

                $rIleGAT += $line[6];
        }
	if ($line[12] =~ /Ile\(TAT\)/){

                $rIleTAT += $line[6];
        }
	if ($line[12] =~ /Leu\(AAG\)/){

                $rLeuAAG += $line[6];
        }
	if ($line[12] =~ /Leu\(CAA\)/){

                $rLeuCAA += $line[6];
        }
	if ($line[12] =~ /Leu\(CAG\)/){

                $rLeuCAG += $line[6];
        }
	if ($line[12] =~ /Leu\(TAA\)/){

                $rLeuTAA += $line[6];
        }
	if ($line[12] =~ /Leu\(TAG\)/){

                $rLeuTAG += $line[6];
        }
	if ($line[12] =~ /Lys\(CTT\)/){

                $rLysCTT += $line[6];
        }
	if ($line[12] =~ /Lys\(TTT\)/){

                $rLysTTT += $line[6];
        }
	if ($line[12] =~ /Met\(CAT\)/){

                $rMetCAT += $line[6];
        }
	if ($line[12] =~ /Phe\(GAA\)/){

                $rPheGAA += $line[6];
        }
	if ($line[12] =~ /Pro\(AGG\)/){

                $rProAGG += $line[6];
        }
	if ($line[12] =~ /Pro\(CGG\)/){

                $rProCGG += $line[6];
        }
	if ($line[12] =~ /Pro\(TGG\)/){

                $rProTGG += $line[6];
        }
	if ($line[12] =~ /Pseudo/){

                $rPseudo += $line[6];
        }
	if ($line[12] =~ /SeC\(e\)\(TCA\)/){

                $rSeC += $line[6];
        }
	if ($line[12] =~ /Ser\(AGA\)/){

                $rSerAGA += $line[6];
        }
	if ($line[12] =~ /Ser\(CGA\)/){

                $rSerCGA += $line[6];
        }
	if ($line[12] =~ /Ser\(TGA\)/){

                $rSerTGA += $line[6];
        }
	if ($line[12] =~ /Ser\(GCT\)/){

                $rSerGCT += $line[6];
        }
	if ($line[12] =~ /Sup/){

                $rSup += $line[6];
        }
	if ($line[12] =~ /Thr\(AGT\)/){

                $rThrAGT += $line[6];
        }
	if ($line[12] =~ /Thr\(CGT\)/){

                $rThrCGT += $line[6];
        }
	if ($line[12] =~ /Thr\(TGT\)/){

                $rThrTGT += $line[6];
        }
	if ($line[12] =~ /Trp\(CCA\)/){

                $rTrpCCA += $line[6];
        }
	if ($line[12] =~ /Tyr\(GTA\)/){

                $rTyrGTA += $line[6];
        }
	if ($line[12] =~ /Undet/){

                $rundet += $line[6];
        }
	if ($line[12] =~ /Val\(AAC\)/){

                $rValAAC += $line[6];
        }
	if ($line[12] =~ /Val\(CAC\)/){

                $rValCAC += $line[6];
        }
	if ($line[12] =~ /Val\(TAC\)/){

                $rValTAC += $line[6];
        }
}
}
print "tRNA\t$tRNA\n";
print "tRNA_reads\t$rtRNA\n";
print "AlaAGC\t$rAlaAGC\n";
print "AlaCGC\t$rAlaCGC\n";
print "AlaTGC\t$rAlaTGC\n";
print "ArgACG\t$rArgACG\n";
print "ArgCCG\t$rArgCCG\n";
print "ArgCCT\t$rArgCCT\n";
print "ArgTCG\t$rArgTCG\n";
print "ArgTCT\t$rArgTCT\n";
print "AsnATT\t$rAsnATT\n";
print "AsnGTT\t$rAsnGTT\n";
print "AspGTC\t$rAspGTC\n";
print "CysGCA\t$rCysGCA\n";
print "GlnCTG\t$rGlnCTG\n";
print "GlnTTG\t$rGlnTTG\n";
print "GluCTC\t$rGluCTC\n";
print "GluTTC\t$rGluTTC\n";
print "GlyCCC\t$rGlyCCC\n";
print "GlyGCC\t$rGlyGCC\n";
print "GlyTCC\t$rGlyTCC\n";
print "HisGTG\t$rHisGTG\n";
print "IleAAT\t$rIleAAT\n";
print "IleGAT\t$rIleGAT\n";
print "IleTAT\t$rIleTAT\n";
print "LeuAAG\t$rLeuAAG\n";
print "LeuCAA\t$rLeuCAA\n";
print "LeuCAG\t$rLeuCAG\n";
print "LeuTAA\t$rLeuTAA\n";
print "LeuTAG\t$rLeuTAG\n";
print "LysCTT\t$rLysCTT\n";
print "LysTTT\t$rLysTTT\n";
print "MetCAT\t$rMetCAT\n";
print "PheGAA\t$rPheGAA\n";
print "ProAGG\t$rProAGG\n";
print "ProCGG\t$rProCGG\n";
print "ProTGG\t$rProTGG\n";
print "SerAGA\t$rSerAGA\n";
print "SerCGA\t$rSerCGA\n";
print "SerTGA\t$rSerTGA\n";
print "SerGCT\t$rSerGCT\n";
print "ThrAGT\t$rThrAGT\n";
print "ThrCGT\t$rThrCGT\n";
print "ThrTGT\t$rThrTGT\n";
print "TrpCCA\t$rTrpCCA\n";
print "TyrGTA\t$rTyrGTA\n";
print "ValAAC\t$rValAAC\n";
print "ValCAC\t$rValCAC\n";
print "ValTAC\t$rValTAC\n";
print "Pseudo\t$rPseudo\n";
print "SeC\t$rSeC\n";
print "SupTTA\t$rSup\n";
print "Undetermined\t$rundet\n";


close $fh;
