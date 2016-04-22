#!/usr/bin/perl
use strict;
use warnings;

open I1, "<", "cytb_accession_number.txt";
open I2, "<", "cytb_sequence.fasta";
open O, ">", "cytb_rename.fa";

my %name;
while (<I1>) {
	chomp;
	my @a = split /\t/;
	$name{$a[1]} = "$a[1]_$a[2]";
}
close I1;

my %fa;
while (<I2>) {
	$/ = ">";
	my $seq = <I2>;
	chomp $seq;
	$seq =~ s/\s//g;
	$/ = "\n";
	/([ABGU]{2}[0-9]{6})/;
	print O ">$name{$1}\n$seq\n";
}
close I2;
close O;

