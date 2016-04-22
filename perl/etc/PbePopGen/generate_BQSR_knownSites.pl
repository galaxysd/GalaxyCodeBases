#!/usr/bin/perl

use strict;
use warnings;

die "This program generate knownSites.vcf for BQSR using reference and BSNP output.\nAuther: Woody
Usage: $0 <reference.fa> <in.nSNP> <in.SNP.gz> [in2.nSNP in2.SNP.gz [...]]\n" if @ARGV < 2;

my $ref_file = shift;
open my $REF, "<", $ref_file;
my %ref;
my %len;
while (<$REF>) {
	chomp;
	s/>//;
	my $seq = <$REF>;
	chomp $seq;
	$ref{$_} = $seq;
	$len{$_} = length $seq;
}
close $REF;

my %v;
foreach (keys %ref) {
	$v{$_} = "0" x $len{$_};
}

while (@ARGV) {
	my %c; #complementary
	foreach (keys %ref) {
		$c{$_} = "n" x $len{$_};
	}
	my $input = shift;
	open my $in, "<", $input;
	warn "Reading $input...\n";
	<$in>;
	while (<$in>) {
		my @a = split / +/;
		substr $c{$a[0]}, $a[1], $a[2], "0" x $a[2];
	}
	close $in;
	$input = shift;
	warn "Reading $input...\n";
	open $in, "-|", "zcat $input";
	<$in>;
	while (<$in>) {
		my @a = split /\t/;
		if ($a[2] eq $a[4]) {
			substr $c{$a[0]}, $a[1], 1, "0";
		}
	}
	close $in;
	foreach my $chr (keys %ref) {
		foreach my $c0 (0 .. $len{$chr}-1) {
			my $base = substr $c{$chr}, $c0, 1;
			if ($base) {
				my $ref_base = substr($ref{$chr}, $c0, 1);
				substr $v{$chr}, $c0, 1, $ref_base unless $ref_base eq "N";
			}
		}
	}
}

warn "Printing ...\n";
print "##fileformat=VCFv4.2\n##reference=file://$ref_file\n";
print "##contig=<ID=$_,length=$len{$_}>\n" foreach sort keys %ref;
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

foreach my $chr (sort keys %ref) {
	foreach my $c0 (0 .. $len{$chr}-1) {
		my $base = substr $v{$chr}, $c0, 1;
		if ($base) {
			my $c1 = $c0 + 1;
			print "$chr\t$c1\t.\t$base\tN\t.\t.\t.\n";
		}	
	}
}
