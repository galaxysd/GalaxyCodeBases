#!/usr/bin/perl
use strict;
use warnings;

open I1, "<", "../a.information/Sample_Info.txt";
my %sam_inf;
while (<I1>) {
	s/^(\S+)//;
	$sam_inf{$1} = $_;
}
close I1;

open I2, "<", "../b.find_hybrid/sample_genotype_region.txt";
open O, ">", "PBE_sample_table.txt";
while (<I2>) {
	chomp;
	/^(\S+)/;
	if ($sam_inf{$1}) {
		print O $_, $sam_inf{$1};
	} else {
		print O $_, "\n";
		warn $1;
	}
}
close I2;
close O;
