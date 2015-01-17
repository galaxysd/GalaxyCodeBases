#!/usr/bin/perl -w
#
# AUthor: Ruan Jue
#
use strict;

my $inf = shift or die("Usage: $0 <ps_fasta> [length_var:20] [min_pe:20] [max_offset:20]\n");
my $var = shift || 20;
my $npe = shift || 20;
my $off = shift || 20;

open(IN, $inf) or die;
my $f = 1;
while(<IN>){
	if(/^>([^_]+)_(\d+)_(\d+)_(\d+)\s/){
		$f = (abs($3 - ($4 % 10000)) <= $var and int(($4 % (100 * 10000)) / 10000) >= $npe and int($4 / (1000 * 10000)) <= $off);
	}
	print if($f);
}
close IN;

1;
