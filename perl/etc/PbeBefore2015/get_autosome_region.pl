#!/usr/bin/perl
use strict;
use warnings;

open I1, "<", "/share/users/miaolin/5.Leopard_cat/3.PBE_paper_2015/c.autosome_individual_structure/a.genotype_from_PBE_paper_2014/SamplesWithAtLeast11loci.txt";
my @sample;
while (<I1>) {
	next if /^Animal/;
	my @a = split /\t/;
	push @sample, $a[0];
}
close I1;

open I2, "<", "/share/users/miaolin/5.Leopard_cat/3.PBE_paper_2015/c.autosome_individual_structure/b.Structure/output_noadm_indep_f/output_noadm_indep_K3_2_f";

while (<I2>) {
	last if /\(%Miss\) :  Inferred clusters/;
}
my @q;
my @region;
while (<I2>) {
	chomp;
	last unless $_;
	my @a = split /:  /;
	$a[1] =~ s/ $//;
	push @q, $a[1];
	my @b = split / /, $a[1];
	if ($b[0] > $b[1] and $b[0] > $b[2]) {
		push @region, "China/Indochina";
	} elsif ($b[1] > $b[0] and $b[1] > $b[2]) {
		push @region, "Southern";
	} elsif ($b[2] > $b[0] and $b[2] > $b[1]) {
		push @region, "Amur";
	} else {
		warn $_;
	}
}
close I2;

die if @sample != @region;

open O, ">", "Sample_AutoRegion.txt";
foreach (0 .. @sample-1) {
	print O "$sample[$_]\t$q[$_]\t$region[$_]\n";
}
close O;
