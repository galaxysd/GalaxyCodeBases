#!/usr/bin/perl
use strict;
use warnings;

my $i = shift;
my $o = shift;

open my $i1, "<", "$i";
open O, ">", "$o";

my $nchar;
my %dna = %{&readfasta($i1)};
close $i1;
my $ntax = keys %dna;

print O "#NEXUS
BEGIN DATA;
Dimensions ntax=$ntax nchar=$nchar;
Format datatype=DNA gap=- missing=?;
Matrix\n";

foreach (sort keys %dna) {
	print O "$_\t$dna{$_}\n";
}
print O ";\nEND;\n";
close O;

sub readfasta {
	my $in = $_[0];
	my %fa;
	while (<$in>) {
		$/ = ">";
		my $seq = <$in>;
		chomp $seq;
		$seq =~ s/\s//g;
		$/ = "\n";
		$nchar = length $seq unless $nchar;
		s/>//;
		/^(\S+)\s/;
		$fa{$1} = $seq;
	}
	return \%fa;
}
