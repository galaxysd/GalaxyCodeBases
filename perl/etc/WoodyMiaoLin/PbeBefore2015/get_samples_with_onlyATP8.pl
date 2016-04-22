#!/usr/bin/perl
use strict;
use warnings;

open I1, "<", "../../../a.mt_haplotype_full_length_tree/b.modified_alignment/haplotype_sampleID.txt";
open I2, "<", "../../a.fromZhangYue/AllPbeATP8_20120315.fasta";
open O, ">", "samples_with_onlyATP8.fasta";

my @s; #sampleIDs with full mt sequence.
while (<I1>) {
	/\t(.+)\n/;
	my @a = split / /, $1;
	push @s, @a;
}
close I1;

while (<I2>) {
	$/ = ">";
	my $seq = <I2>;
	chomp $seq;
	$seq =~ s/\s//g;
	$/ = "\n";
	s/>//;
	/^(\S+)\s/;
	my $b = $1;
	if (!grep(/^$b$/, @s)) {
		print O ">$1\n$seq\n";
	}
}
close I2;
close O;
