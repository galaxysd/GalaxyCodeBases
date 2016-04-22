#!/usr/bin/perl
use strict;
use warnings;

open I1, "<", "sample_location.txt";
my %sl;
while (<I1>) {
	/(\w+)\t(\w+)\n/;
	$sl{$1} = $2;
}
close I1;

open I2, "<", "../b.haplotype/PbeATP8_haplotype.txt";
open O, ">", "PbeATP8_haplotype_location_sample.txt";
while(<I2>) {
	chomp;
	my @a = split /\t/;
	my @b = split / /, "$a[1]";
	my %lo;
	foreach (@b) {
		if ($sl{$_}) {
			$lo{$sl{$_}} = 1;
			$sl{$_} = 0;
		} else {
			warn $_;
		}
	}
	print O "$a[0]\t";
	print O join(" ", sort keys %lo);
	print O "\t$a[1]\n";
	my $c = <I2>;
	print O $c;
}
close I2;
close O;

foreach (sort keys %sl) {
	warn $_ if $sl{$_};
}
