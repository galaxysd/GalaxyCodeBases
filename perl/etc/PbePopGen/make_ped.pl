#!/usr/bin/perl
use strict;
use warnings;

die "This program makes PED file using gzipped fasta files.\nAuther: Woody
Usage: $0 <list (column1: sample, column2: file, column3: 1=male|2=female, column4: population)> <output prefix>\n" if @ARGV < 2;

open LIST, "<", "$ARGV[0]";
open GT, ">", "$ARGV[1].ped";
open SNP, ">", "$ARGV[1].pedsnp";
open IND, ">", "$ARGV[1].pedind";

my @id;
my %sex;
my %pop;
my %file;
while (<LIST>) {
	chomp;
	my @a = split / /;
	push @id, $a[0];
	$sex{$a[0]} = $a[2];
	$pop{$a[0]} = $a[3];
	open $file{$a[0]}, "-|", "zcat $a[1]";
}
close LIST;




my %chr_gt;
my %chr_coo;
while (readline $file{$id[0]}) {
	chomp;
	s/^>chr//;
	last unless /^\w\d$/;
	my $chr0 = $_;
	my $seq0 = readline $file{$id[0]};
	chomp $seq0;
	my $len = length $seq0;
	my @seq;
	push @seq, $seq0;
	for (@id) {
		next if $_ eq $id[0];
		my $chr1 = readline $file{$_};
		chomp $chr1;
		warn if $chr1 ne ">chr$chr0";
		my $seq1 = readline $file{$_};
		push @seq, $seq1;
	}
	for my $coo (0 .. $len-1) {
		my @base;
		my %base;
		for (@seq) {
			my $b = substr($_, $coo, 1);
			push @base, $b;
			++$base{$b};
		}
		my $nb = keys %base;
		if ($base{N}) {
			next;
#			next if $nb == 2;
#			next if $base{N} > 10;
		} else {
			next unless $nb == 2;
		}
		for (0 .. @id-1) {
			$chr_gt{$chr0}[$_] .= $base[$_];
		}
		push @{$chr_coo{$chr0}}, $coo;
	}
}

my %base_gt = ("A","A A", "C","C C", "G","G G", "T","T T", "R","A G", "Y","C T", "M","A C", "K","G T", "S","G C", "W","A T", "N","0 0");
for my $i (0 .. @id-1) {
	close $file{$id[$i]};
	print IND "0\t$id[$i]\t0\t0\t$sex{$id[$i]}\t$pop{$id[$i]}\n";
	print GT "0\t$id[$i]\t0\t0\t$sex{$id[$i]}\t$pop{$id[$i]}";
	for my $c (sort keys %chr_gt) {
		my $len = length $chr_gt{$c}[$i];
		for (0 .. $len-1) {
			my $b = substr $chr_gt{$c}[$i], $_, 1;
			print GT "\t$base_gt{$b}";
		}
	}
	print GT "\n";
}
close IND;
close GT;
my %chr_num = qw / A1 1 A2 2 A3 3 B1 4 B2 5 B3 6 B4 7 C1 8 C2 9 D1 10 D2 11 D3 12 D4 13 E1 14 E2 15 E3 16 F1 17 F2 18 /;
for my $c (sort keys %chr_coo) {
	for (0 .. @{$chr_coo{$c}}-1) {
		print SNP "$chr_num{$c}\tchr${c}snp";
		printf SNP "%08d", 1+$_;
		print SNP "\t0\t$chr_coo{$c}[$_]\n";
	}
}
close SNP;
