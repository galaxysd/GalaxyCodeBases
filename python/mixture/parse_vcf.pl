#!/usr/bin/env perl
use strict;
use warnings;

use Data::Dump qw(ddx);

my $Usage = "Usage: $0 <vcf.gz>\n";
die $Usage if @ARGV < 1;

my ($filename) = @ARGV;
#my $cmd = "bcftools view -m3 -v snps $filename | bcftools query -f '%CHROM\t%POS\t%REF,%ALT\t%QUAL[\t%TGT:%AD:%GQ]\n'";
my $cmd = "bcftools view -U -m2 -i '%QUAL>=40 & MIN(FMT/GQ)>20 & GT[2]=\"het\"' -v snps $filename | bcftools query -f '%CHROM\t%POS\t%REF,%ALT\t%QUAL[\t%TGT:%AD:%GQ]\n'";
=pod
chr1	102951	C,T	1575.58	C/C/C/T:23,6:25	C/C/C/T:88,27:67	C/C/C/T:111,33:93
chr1	633963	C,T	6936.96	T/T/T/T:5,242:99	C/C/C/C:250,0:99	C/C/C/T:172,77:75
chr1	889636	A,G,T	1356.65	A/A/G/T:14,3,8:21	A/A/G/T:36,20,10:48	A/A/G/T:50,23,18:75
chr1	1099396	G,A,C	214.66	G/G/G/G:22,1,0:24	G/G/A/C:67,5,9:58	G/G/A/C:89,6,9:24
=cut
my @SampleIDs;
open(SID,'-|',"bcftools query -l $filename") or die "Error opening [$filename]: $!\n";
while(<SID>) {
	chomp;
	push @SampleIDs,$_;
}
close SID;
ddx \@SampleIDs;
# ["T10C", "T3C", "mixed"] => v,k,m
open(IN,"-|",$cmd) or die "Error opening [$filename]: $!\n";
while (<IN>) {
	print $_;
	chomp;
	my ($Chrom,$Pos,$RefAlt,$Qual,@GTs) = split /\t/,$_;
	my %GTs;
	for my $i (0 .. $#SampleIDs) {
		my ($sGT,$sAD,$sGQ) = split /:/,$GTs[$i];
		my @aGT = split /[|\/]/,$sGT;
		my @aAD = split /\,/,$sAD;
		my %counter;
		for (@aGT) {
			++$counter{$_};
		}
		#my @result = sort(keys %counter);
		$GTs{$SampleIDs[$i]} = \%counter;
	}
	$GTs{'_R'} = {%{$GTs{'mixed'}}};
	for (keys %{$GTs{'T10C'}}) {
		if (exists $GTs{'_R'}{$_}) {
			delete $GTs{'_R'}{$_};
		}
	}
	ddx \%GTs;
}
close IN;
