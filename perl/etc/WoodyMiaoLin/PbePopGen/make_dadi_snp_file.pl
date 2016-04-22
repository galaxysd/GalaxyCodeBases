#!/usr/bin/perl
use strict;
use warnings;

die "This program makes DaDi SNP data format file using gzipped fasta files.\nAuther: Woody
Usage: $0 <list (column1: sample, column2: file, column3: \"Column1\"|\"Column2\"|population)> <output>\n" if @ARGV < 2;

open LIST, "<", $ARGV[0];
open O, ">", $ARGV[1];

my $C1;
my $C2;
my @id;
my @po;
my %pop;
my %id_pop;
my %file;
while (<LIST>) {
	chomp;
	my @a = split / /;
	open $file{$a[0]}, "-|", "zcat $a[1]";
	if ($a[2] eq "Column1") {
		$C1 = $a[0];
	} elsif ($a[2] eq "Column2") {
		$C2 = $a[0];
	} else {
		push @id, $a[0];
		push @po, $a[2];
		$pop{$a[2]} = 1;
	}
}
close LIST;

my %base_allele = qw / A AA  C CC  G GG  T TT  R AG  Y CT  M AC  K GT  S GC  W AT  N NN /;
my %chr_gt;
print O "$C1\t$C2\tAllele1\t", join("\t", sort keys %pop), "\tAllele2\t", join("\t", sort keys %pop), "\tChromosome\tOneBasedCoordinate\n";
while (readline $file{$C1}) {
	chomp;
	s/^>//;
	last unless /^chr\w\d$/;
	my $chr0 = $_;
	my $seq0 = readline $file{$C1};
	chomp $seq0;
	my $len = length $seq0;
	my @seq;
	push @seq, $seq0;
	readline $file{$C2};
	$seq0 = readline $file{$C2};
	chomp $seq0;
	push @seq, $seq0;
	for (@id) {
		my $chr1 = readline $file{$_};
		chomp $chr1;
		warn if $chr1 ne ">$chr0";
		my $seq1 = readline $file{$_};
		push @seq, $seq1;
	}
	for my $coo (0 .. $len-1) {
		my @base;
		my %base;
		for (2 .. @seq-1) {
			my @b = split //, $base_allele{substr($seq[$_], $coo, 1)};
			push @base, \@b;
			++$base{$b[0]};
			++$base{$b[1]};
		}
		if ($base{N}) {
			next;
		} else {
			next if keys %base != 2;
		}
		my ($allele0, $allele1) = keys %base;
		my %pop_allele;
		for (0 .. @base-1) {
			++$pop_allele{$po[$_]}{$base[$_][0]};
			++$pop_allele{$po[$_]}{$base[$_][1]};
		}
		my @num_allele0;
		my @num_allele1;
		for (sort keys %pop) {
			$pop_allele{$_}{$allele0} = 0 unless $pop_allele{$_}{$allele0};
			$pop_allele{$_}{$allele1} = 0 unless $pop_allele{$_}{$allele1};
			push @num_allele0, $pop_allele{$_}{$allele0};
			push @num_allele1, $pop_allele{$_}{$allele1};
		}
		print O substr($seq[0], $coo-1, 3), "\t", substr($seq[1], $coo-1, 3), "\t$allele0\t";
		print O join("\t", @num_allele0);
		print O "\t$allele1\t";
		print O join("\t", @num_allele1);
		print O "\t$chr0\t", $coo+1, "\n";
	}
}
close O;
