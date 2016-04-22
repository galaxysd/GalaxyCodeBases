#!/usr/bin/perl
use strict;

open I, "<", "a2.txt";
open O, ">", "a3.txt";

my $head = <I>;
$head =~ s/\r//;
print O $head;

while (<I>) {
	s/\r\n//;
	my @a = split /\t/;
	my $group = $a[5].$a[8].$a[11];
	if (@{[$group=~/North/g]} >= 2) {
		if ($group =~ /Sunda/) {
			$a[12] = "n/a";
		} else {
			$a[12] = "North";
		}
	} elsif (@{[$group=~/Sunda/g]} >= 2) {
		if ($group =~ /North/) {
			$a[12] = "n/a";
		} else {
			$a[12] = "Sunda";
		}
	} else {
		$a[12] = "n/a";
	}
	print O join("\t", @a), "\n";
}
close I;
close O;
