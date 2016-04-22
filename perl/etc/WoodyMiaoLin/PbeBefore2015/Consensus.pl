#!/usr/bin/perl
use strict;
use warnings;

die("
Program: Consensus
Version: 1.0
Release: Aug. 26, 2014\n
Auther: Woody
Usage: $0 <input> <output>\n
Both input and output files are fasta format.
Make consensus from Squencher \"File>Export>Contig\". Partial sequences must have identical ID.
\n") if (@ARGV<2);

my $in = shift;
my $out = shift;
open I, "<", "$in";

# Read input file.
my %input; # input sequences
my $nchar; # number of characters
my $nin; # number of input sequences
while (<I>) {
	$/ = ">";
	my $seq = <I>;
	chomp $seq;
	$seq =~ s/\s//g;
	$/ = "\n";
	$nchar = length $seq unless $nchar;
	s/>//;
	/^(\S+)\s/;
	push @{$input{$1}}, $seq;
	++$nin;
}
close I;
warn "Read input complete!\n";
warn "No. sequences = $nin\n";

# Make sample consensus sequence.
my %samseq; # sample consensus sequences
my @samid; # sample IDs
foreach my $a (sort keys %input) {
	push @samid, $a;
	my $con = $input{$a}[0];
	my $n = @{$input{$a}} - 1;
	if ($n == 0) {
		$samseq{$a} = $con;
	} else {
		foreach my $b (1 .. $n) {
			for (my $i = 0; $i < $nchar; $i++) {
				my $j = $i + 1;
				if (substr($input{$a}[$b],$i,1) eq "-") {
					next;
				} else {
					if (substr($con,$i,1) eq "-") {
						substr($con,$i,1) = substr($input{$a}[$b],$i,1);
					} else {
						if (substr($con,$i,1) eq substr($input{$a}[$b],$i,1)) {
							next;
						} else {
							my @nt = (substr($con,$i,1), substr($input{$a}[$b],$i,1));
							if (grep(/A/, @nt) and grep(/G/, @nt)) {
								substr($con,$i,1) = "R";
							} elsif (grep(/C/, @nt) and grep(/T/, @nt)) {
								substr($con,$i,1) = "Y";
							} elsif (grep(/A/, @nt) and grep(/C/, @nt)) {
								substr($con,$i,1) = "M";
							} elsif (grep(/G/, @nt) and grep(/T/, @nt)) {
								substr($con,$i,1) = "K";
							} elsif (grep(/C/, @nt) and grep(/G/, @nt)) {
								substr($con,$i,1) = "S";
							} elsif (grep(/A/, @nt) and grep(/T/, @nt)) {
								substr($con,$i,1) = "W";
							} else {
								next;
							}
						}
					}
				}
			}
		}
		$samseq{$a} = $con;
	}
}
open O1, ">", $out;
print O1 ">$_\n$samseq{$_}\n" foreach sort keys %samseq;
close O1;
warn "Make sample consensus complete and print to $out!\n";
my $aa1 = keys %samseq;
warn "No. samples = $aa1\n";
