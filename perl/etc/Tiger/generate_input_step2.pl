#!/usr/bin/perl
use strict;
use warnings;

open I, "<", "scaffold75_1458SNP.16.gt";
open O, ">", "scaffold75_1458.16.inp";

my $id = <I>;
chomp $id;
my @id = split /\t/, $id;
my (@pos, @gt);
while (<I>) {
	chomp;
	my @a = split /\t/;
	push @pos, $a[1];
	foreach (2 .. 17) {
		my @b = split /\//, $a[$_];
		$gt[$_][0] .= $b[0];
		$gt[$_][1] .= $b[1];
	}
}
close I;

my $num_indi = @id-2;
my $num_loci = @pos;
print O "$num_indi\n";
print O "$num_loci\n";
print O (join " ", "P", @pos), "\n";
print O ("S" x $num_loci), "\n";
foreach (2 .. 17) {
	print O $id[$_], "\n";
	print O $gt[$_][0], "\n";
	print O $gt[$_][1], "\n";
}
close O;
