#!/usr/bin/perl
use strict;
use warnings;

open O, ">", "2.2.psmc.makefile";

my @a = ("PBE084", "PBE144", "PVI033");
my @target;
foreach my $i (1 .. 100) {
	foreach my $a (@a) {
		push @target, "$a.round-$i.psmc";
		print O "$a.round-$i.psmc:\n";
		print O "\tpsmc -r 5 -p \"4+25*2+4+6\" -b $a.split.psmcfa >$a.round-$i.psmc\n"
	}
}
print O "targets: ", join(" ", @target), "\n";
close O;
