#!/usr/bin/perl
use strict;
use warnings;

my $i = shift;
my $o = shift;

open I, "<", "$i";
open O, ">", "$o";

while (<I>) {
	last if /matrix/i;
}

while (<I>) {
	last if /;/;
	my @a = split /\t/;
	print O ">$a[0]\n$a[1]";
}

close I;
close O;
