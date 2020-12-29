#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dump qw(ddx);

my $s1HomoLevel = 0.8;
my $s2HomoLevel = 0.3;
my $s2HeteLevel = 0.2;
mu $minAllelDepth = 100;

my %StrDat;
while(<>) {
	chomp;
	my ($aID,$depth) = split /\s+/;
	my ($strID,$strAllele) = split /\-/,$aID;
	#print join("\t",$strID,$strAllele,$depth),"\n";
	$StrDat{$strID}{$strAllele} = $depth;
}

ddx \%StrDat;
