#!/bin/env perl
use strict;
use warnings;

unless (@ARGV){
	print "perl $0 <in>\n";
	exit;
}

my ($fq1) = @ARGV;
open FQ,'-|',"gzip -dc $fq1";

while (<FQ>) {
	my $a=$_;
	my $b=<FQ>;
	my $c=<FQ>;
	my $d=<FQ>;
	$b=substr $b,0,50;
	$d=substr $d,0,50;
	print "$a$b\n$c$d\n";
}
