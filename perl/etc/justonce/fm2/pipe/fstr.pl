#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dump qw(ddx);

my $s1HomoLevel = 0.8;
my $s2HomoLevela = 0.75;
my $s2HeteLevela = 0.3;
my $s2HeteLevelb = 0.2;
my $minAllelDepth = 100;

my %StrDat;
while(<>) {
	chomp;
	my ($aID,$depth) = split /\s+/;
	my ($strID,$strAllele) = split /\-/,$aID;
	#print join("\t",$strID,$strAllele,$depth),"\n";
	$StrDat{$strID}{$strAllele} = $depth;
}
#ddx \%StrDat;

sub SumArryItr { my $agg = 0; $agg += $_ for @_;  return $agg }

for my $k (keys %StrDat) {
	my %oneSTR = %{$StrDat{$k}};
	#ddx \%oneSTR;
	my @Allels = sort {$oneSTR{$b} <=> $oneSTR{$a}} keys %oneSTR;
	#ddx \@Allels;
	my $Sum = SumArryItr(values %oneSTR);
	#print $Sum,"\n---\n";
	my @GT;
	if ($oneSTR{$Allels[0]} <= 100) {
		@GT = ('NA');
	} elsif ($oneSTR{$Allels[0]} >= $Sum * $s1HomoLevel ) {
		@GT = ($Allels[0]);
	} elsif ( ($oneSTR{$Allels[0]} >= $Sum * $s2HomoLevela) and ($oneSTR{$Allels[1]} < $Sum * $s2HeteLevelb) ) {
		@GT = ($Allels[0]);
	} elsif ( ($oneSTR{$Allels[0]} >= $Sum * $s2HeteLevela) and ($oneSTR{$Allels[1]} >= $Sum * $s2HeteLevelb) ) {
		@GT = ($Allels[0],$Allels[1]);
	} else {
		@GT = ('Unknown');
	}
	ddx \@GT;
}



