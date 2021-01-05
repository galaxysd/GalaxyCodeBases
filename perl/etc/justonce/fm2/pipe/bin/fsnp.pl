#!/usr/bin/env perl

use strict;
use warnings;

#use Data::Dump qw(ddx);

my $minQual = 20;
my $minDepth = 50;

my %SnpDat;
while(<>) {
	chomp;
	my ($rsID,$depth,$qual,$TGT) = split /\s+/;
	next if $rsID !~ '^rs';
	#print join("\t",$rsID,$depth,$qual,$TGT),"\n";
	if ($depth <= $minDepth or $qual <= $minQual) {
		$TGT = '_NA_';
	}
	print "SNP\t$rsID: $TGT\n";
}
#ddx \%StrDat;
