#!/usr/bin/perl
use strict;
use warnings;

open I1, "<", "Sample_Info.txt";
open I2, "<", "Origin_STR_SNP.txt";
open O, ">", "STR_gt_info.txt";

my %info;
while (<I1>) {
	s/\r\n$//;
	my @a = split /\t/;
	$info{$a[0]} = $_;
}
close I1;

while (<I2>) {
	s/\r\n$//;
	my @a = split /\t/;
	print O $_, "\t", $info{$a[0]}, "\n";
}
close I1;
close I2;
close O;
