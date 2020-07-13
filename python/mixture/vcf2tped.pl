#!/usr/bin/env perl
use strict;
use warnings;

use Data::Dump qw(ddx);

my $Usage = "Usage: $0 <vcf.gz> <out prefix>\n";
die $Usage if @ARGV < 1;

my ($filename,$outp) = @ARGV;
my $cmd = "bcftools view -U -m2 -i '%QUAL>=40 & MIN(FMT/GQ)>20 & GT[2]=\"het\"' -v snps $filename |bcftools view -e 'FMT/DP=\".\"'| bcftools query -f '%CHROM\t%POS\t%REF,%ALT\t%QUAL[\t%TGT:%AD:%GQ]\n'";
my @SampleIDs;
open(SID,'-|',"bcftools query -l $filename") or die "Error opening [$filename]: $!\n";
while(<SID>) {
	chomp;
	push @SampleIDs,$_;
}
close SID;
ddx \@SampleIDs;

open O,'>',"$outp.tfam" or die "[$outp.tfam]:$!\n";
for (@SampleIDs) {
	my $FID = 'T';
	my $IID = $_;
	print O join(' ',$FID,$IID,0,0,0,0),"\n";
	# https://www.cog-genomics.org/plink/1.9/formats#fam
}
close O;

open O,'>',"$outp.tped" or die "[$outp.tped]:$!\n";

# https://www.cog-genomics.org/plink/1.9/formats#tped
close O;
