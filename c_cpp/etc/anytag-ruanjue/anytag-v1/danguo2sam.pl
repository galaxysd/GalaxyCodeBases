#!/usr/bin/perl -w
#
#Author: Ruan JUe
#
use strict;

my $seqf = shift or die("Usage: $0 <fasta_file> <danguo_file.gz>\n");
my $danf = shift or die("Usage: $0 <fasta_file> <danguo_file.gz>\n");

open(IN, $seqf) or die;
my $sid;
my $sln;
while(<IN>){
	if(/^>(\S+)/){
		print "\@SQ SN:$sid LN:$sln\n" if($sid);
		$sid = $1;
		$sln = 0;
	} else {
		chomp;
		$sln += length $_;
	}
}
print "\@SQ SN:$sid LN:$sln\n" if($sid);
close IN;

open(IN, "gzip -dc $danf |") or die;
my $ref;
my $default_q = chr(ord('B') + 20);
while(<IN>){
	if(/^B (\d+)/){
		$ref = "longread$1";
	} elsif(/^I (\d+)\t([+-])\t(\S+)\t(\S+)/){
		my $rid = "R$1";
		my $flag = ($2 eq '-')? 0x10 : 0x0;
		my $cigar = $3;
		my $seq = $4;
		my $qual = $default_q x length($seq);
		my $pos = 0;
		if($cigar=~/^(\d+)D/){
			$pos = $1;
		}
		$cigar=~s/^\d+D//;
		$cigar=~s/\d+D$//;
		print "$rid\t$flag\t$ref\t$pos\t0\t$cigar\t$seq\t$qual\tXT:A:U\n";
	}
}
close IN;

1;
