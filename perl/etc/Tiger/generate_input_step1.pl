#!/usr/bin/perl
use strict;
use warnings;

open I1, "<", "scaffold75SNP.18.gt";
open I2, "<", "scaffold1458SNP.18.gt";
open O, ">", "scaffold75_1458SNP.16.gt";

print O "#Scaffold\tPosition\tJHH001\tGZXJ03\tGZXJ05\tGZXJ26\tGZXJ27\tGZXJ28\tGZXJ29\tGZXJ30\tGZXJ33\tBHX011\tBHX019\tGZXJ04\tGZXJ06\tGZXJ31\tGZXJ32\tGZXJ38\n";

while (<I1>) {
	chomp;
	my @a = split /\t/;
	foreach (4 .. 21) {
		my @b = split /\//, $a[$_];
		$a[$_] = "$a[2]/$a[2]" if (($b[0] eq 0) and ($b[1] eq 0));
		$a[$_] = "$a[3]/$a[3]" if (($b[0] eq 1) and ($b[1] eq 1));
		$a[$_] = "$a[2]/$a[3]" if (($b[0] eq 0) and ($b[1] eq 1));
		$a[$_] = "$a[3]/$a[2]" if (($b[0] eq 1) and ($b[1] eq 0));
	}
	splice @a, -3, 2;
	splice @a, 2, 2;
	print O (join "\t", @a), "\n";
}

while (<I2>) {
       	chomp;
	my @a = split /\t/;
	$a[1] += 5658748;
        foreach (4 .. 21) {
		my @b = split /\//, $a[$_];
		$a[$_] = "$a[2]/$a[2]" if (($b[0] eq 0) and ($b[1] eq 0));
		$a[$_] = "$a[3]/$a[3]" if (($b[0] eq 1) and ($b[1] eq 1));
		$a[$_] = "$a[2]/$a[3]" if (($b[0] eq 0) and ($b[1] eq 1));
		$a[$_] = "$a[3]/$a[2]" if (($b[0] eq 1) and ($b[1] eq 0));
	}
	splice @a, -3, 2;
	splice @a, 2, 2;
	print O (join "\t", @a), "\n";
}

close I1;
close I2;
close O;
