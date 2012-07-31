#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Data::Dump qw(ddx);
use Galaxy::IO;
use Galaxy::IO::FASTAQ qw(readfq getQvaluesFQ);
use Galaxy::SeqTools;

die "Usage: $0 <ec list> <bcgv bcf> <out>\n" if @ARGV<2;
# no more "<sorted bam files>", use `samtools merge -r <out.bam> <in1.bam> <in2.bam> [...]` first !
my $eclst=shift;
my $bcfs=shift;
my $outfs=shift;

my $ReadLen = 101-5;
my $minOKalnLen = 30;
my $minAlignLen = int($ReadLen * 0.6);
$minAlignLen = $minOKalnLen if $minAlignLen < $minOKalnLen;
die if $ReadLen < $minOKalnLen;

my $Eseq="CTGCAG";
my $EseqLen = length $Eseq;
my $EcutAt=5;
my $EseqL="CTGCA";#substr $Eseq,0,$EcutAt;
my $EseqR="TGCAG";

my $EfwdTerm5=1-$EcutAt+1;	# -3
my $EfwdTerm3=1-$EcutAt+$ReadLen;	# 92
my $ErevTerm3=0;	# 0
my $ErevTerm5=$ErevTerm3-$ReadLen+1;	# -95
# 12008 -> [12005,12099],[11914,12008]
my $Rfwd2Ec = -$EfwdTerm5;	# 3
my $Rrev2Ec = -$ErevTerm3;	# 0

# C|TGCAG
my $PosLeft = -$EfwdTerm3;	# -92
my $PosRight = -$ErevTerm5;	# 95;
my $PosECsft = $Rfwd2Ec;	# 3;
# - [x-$PosLeft,x+$PosECsft], Len=96
# + [x,x+$PosRight], Len=96
# A [x-$PosLeft,x+$PosRight], Len=95+92+1=186

my %Stat;
my $t;
open O,'>',$outfs or die "Error opening $outfs : $!\n";
$t = "# EClst: [$eclst], Enzyme: [$Eseq], Cut after $EcutAt\n# Bams: [$bcfs]\n\n";
print O $t;
print $t;

my %Markers;
open L,'<',$eclst or die;
while (<L>) {
	next if /^(#|$)/;
	my ($chr,$pos,$strand,$mark,$count,$samples) = split /\t/;
#warn "$chr,$pos,$strand,$samples,$mark\n";
	next if $samples < 2;
	push @{$Markers{$chr}},[$pos,$strand]
}
close L;

my $th = openpipe('bcftools view',$bcfs);
my @Samples;
while (<$th>) {
	next if /^##/;
	my @data = split /\t/;
	if ($data[0] eq '#CHROM') {
		@Samples = map {my $t=(split /\//)[-1];$t=~s/_cut//g;$t=~s/-/./g; $_=join('_',(split /\./,$t)[-5,-6]);} splice @data,9;
		# ../5.bam_0000210210_merged/d1_4_merged.D4.JHH001.XTU.sort.rmdup.bam
		last;
	}
}

__END__
my %VCF;
while (<$th>) {
	next if /^#/;
	my @data = split /\t/;
}