#!/usr/bin/perl
use strict;
use warnings;

open I, "-|", "gzip -dc ../4.G-PhoCS/input.gz";
open O, "|-", "gzip -9c >GPhoCS_PbePvi_input.gz";

while (<I>) {
	if (/^chr/) {
		my @a = split /\t/;
		print O "$a[0]\t2\t$a[2]";
	} elsif (/^FCA|PVI/) {
		next;
	} else {
		print O $_;
	}
}
close I;
close O;
