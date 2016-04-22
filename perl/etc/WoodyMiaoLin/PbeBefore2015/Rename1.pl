#!/usr/bin/perl
use strict;
use warnings;

my $in = shift;

open I, "<", $in;
open O, ">", "rename_$in";

while (<I>) {
	/((Pbe|Pvi|Ipl|Mir|Ppa|Pti|Pte|Pma)\d+)/;
	my $name = $1;
	$/ = ">";
	my $seq = <I>;
	chomp $seq;
	$seq =~ s/\s//g;
	$/ = "\n";
	warn $_ unless $name;
	print O ">$name\n$seq\n";
}
close I;
close O;
