#!/usr/bin/perl
use strict;
use warnings;

open I, "<", "../b.modified_alignment/haplotype_sampleID.txt";
open O, ">", "tip_label.txt";

my $total;
while (<I>) {
	chomp;
	my @a = split /\t/;
	my @b = split / /, $a[1];
	my $prefix;
	my @number;
	my $n;
	foreach (@b) {
		++$n;
		/(\D+)(\d+)/;
		if ($prefix) {
			warn $_ if $prefix ne $1;
		} else {
			$prefix = $1;
		}
		push @number, $2;
	}
	$total += $n;
	my $c = $prefix . join(",", sort {$a <=> $b} @number) . ", n = $n";
	print O "$a[0]\t($c)\n";
}
warn $total;
close I;
close O;
