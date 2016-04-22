#!/usr/bin/perl
use strict;
use warnings;

my $fn1 = shift;
$fn1 =~ m/^([.-z]+)R1_001.fastq.gz$/;
my $fn2 = "${1}R1_001.fastq.gz";
$fn1 =~ m/(\w+)_R1_001.fastq.gz$/;
my $f = $1; # file prefix;

open I1, "-|", "cutadapt -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATC -g AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -e 0.08 -n 2 -O 5 $fn1 2>$f.1.cut.log"  or die "Error opening $fn1: $!\n";
open I2, "-|", "cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -g GATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT -e 0.08 -n 2 -O 5 $fn2 2>$f.2.cut.log"  or die "Error opening $fn2: $!\n";
open O1, "|-", "gzip -9c >$f_CutAdd_1.fq.gz";
open O2, "|-", "gzip -9c >$f_CutAdd_2.fq.gz";
open LOG, ">>", "count_fq_n_addN.log";

my $n1 = 0; # add N number;
my $b1 = 0; # basepair number;
my $r1 = 0; # reads number;
my $n2 = 0;
my $b2 = 0;
my $r2 = 0;

while (<I1>) {
	my $r11 = $_;
	my $r12 = <I1>;
	my $r13 = <I1>;
	my $r14 = <I1>;
	my $l1 = length $r12;
	if ($l1 == 1) {
		print O1 $r11;
		print O1 "N\n";
		print O1 $r13;
		print O1 "#\n";
		++$n1;
		++$r1;
	} elsif ($l1 > 1) {
		print O1 $r11;
		print O1 $r12;
		print O1 $r13;
		print O1 $r14;
		--$l1;
		$b1 += $l1;
		++$r1;
	}
	my $r21 = <I2>;
	my $r22 = <I2>;
	my $r23 = <I2>;
	my $r24 = <I2>;
	my $l2 = length $r22;
	if ($l2 == 1) {
		print O2 $r21;
		print O2 "N\n";
		print O2 $r23;
		print O2 "#\n";
		++$n2;
		++$r2;
	} elsif ($l2 > 1) {
		print O2 $r21;
		print O2 $r22;
		print O2 $r23;
		print O2 $r24;
		--$l2;
		$b2 += $l2;
		++$r2;
	}
}
close I1;
close I2;
close O1;
close O2;

my $r;
if ($r1 == $r2) {
	$r = $r1 + $r2;
} else {
	$r = "error";
}
my $b = $b1 + $b2;
my $n = $n1 + $n2;

print LOG "$f\t$n1\t$b1\t$r1\t$n2\t$b2\t$r2\t$n\t$b\t$r\n";
close LOG;
