#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 1.0.0 @ 20180417
=cut
use strict;
use warnings;
use Galaxy::IO;
use Data::Dump qw(ddx);

my $minDepth = 100;

#my $in = openfile('MA20.snp.501a100.txt.gz');
my $in = openpipe("bcftools query -f '%CHROM\t%POS\t%REF,%ALT\t%QUAL[\t%TGT;%AD]\n' -i'POS=501 && AVG(FMT/AD)>=$minDepth'",'MA20.snp.gz');
my @Samples = qw(C1 F1 F2 F3 F4 F5 M1 M2 M3 M4 M5);
my (%Result,%GroupTGT);

while(<$in>) {
	chomp;
	my ($chr,$pos,$Alleles,$Qual,@SampleDat) = split /\t/;
	#print join(", ",$chr,$pos,@SampleDat),"\n";
	my (@SampleGT,@SampleDep);
	for (@SampleDat) {
		my ($TGT,$Deps) = split /;/;
		my @aSampleGT = split /[\/|]/,$TGT;
		my @aSampleDep = split /,/,$Deps;
		my $aSumDep;
		map { $aSumDep += $_ } @aSampleDep;
		push @SampleDep,$aSumDep;
		push @SampleGT,join('',sort @aSampleGT);
	}
	my (%GTf,%GTm);
	for (1..5) {
		++$GTf{$SampleGT[$_]} if $SampleDep[$_] >= $minDepth;
	}
	for (6..10) {
		++$GTm{$SampleGT[$_]} if $SampleDep[$_] >= $minDepth;
	}
	#++$GTf{$SampleGT[$_]} for (1..5);
	#++$GTm{$SampleGT[$_]} for (6..10);
	++$Result{'F'}{keys(%GTf)};
	++$Result{'M'}{keys(%GTm)};
	if (keys(%GTf)+keys(%GTm)>2) {
		print "@SampleGT,@SampleDep\n";
		ddx %GTf if keys(%GTf)>1;
		ddx %GTm if keys(%GTm)>1;
	}
	#print join(", ",$chr,$pos,@SampleGT,@SampleDep),"\n";
	#ddx \%Result;
}
close $in;

ddx \%GroupTGT;
ddx \%Result;
