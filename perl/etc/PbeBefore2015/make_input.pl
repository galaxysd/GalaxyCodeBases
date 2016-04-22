#!/usr/bin/perl
use strict;
use warnings;

open I, "<", "Origin_STR_SNP_All.txt";
open O1, ">", "Origin_STR_SNP_All_12loci.txt";
open O2, ">", "StructureInput.txt";

my $a = <I>;
$a =~ s/\r//;
print O1 $a;
while (<I>) {
	s/\r\n//;
	my @a = split /\t/;
	splice @a, 0, 3;
	my $n = 0;
	foreach (0 .. 33) {
		++$n if (int($_/2) == $_/2) and $a[$_];
		if ($_ >= 20 and $_ <= 29) {
			$a[$_] =~ s/^1//;
		}
	}
	if ($n >= 12) {
		print O1 "$_\n";
		print O2 join("\t", @a), "\n";
	}
}
close I;
close O1;
close O2;
