#!/usr/bin/perl
use strict;
use warnings;

die "This program picks up loci with certain ancestral allele from DaDi SNP data format file, and arranges Allele1 to be ancestral.\nAuther: Woody
Usage: $0 <dadi snp file> <output>\n" if @ARGV < 2;

open I, "<", $ARGV[0];
open O, ">", $ARGV[1];

my $head = <I>;
print O $head;
my @h = split /\t/, $head;
my $index_a2;
for (0 .. @h-1) {
	if ($h[$_] eq "Allele2") {
		$index_a2 = $_;
		last;
	}
}
my $index_lastpop = 2 * ($index_a2 - 2) + 1;

while (<I>) {
	my @a = split /\t/;
	my @o1 = split //, $a[0];
	next if $o1[1] eq "N";
	my @o2 = split //, $a[1];
	next if $o1[1] ne $o2[1];
	if ($o1[1] eq $a[2]) {
		print O $_;
	} elsif ($o1[1] eq $a[$index_a2]) {
		print O join "\t", @a[0, 1, $index_a2..$index_lastpop, 2..$index_a2-1, $index_lastpop+1..@a-1];
	} else {
		next;
	}
}
close I;
close O;
