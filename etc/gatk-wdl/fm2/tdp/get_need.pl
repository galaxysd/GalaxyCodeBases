#!/usr/bin/env perl
use strict;

my ($need,$all) = @ARGV;

my %freq;
open ALL,"<$all" or die($!);
while (<ALL>){
	chomp;
	my @data = split /\t/,$_;
	$freq{$data[0]}{RG} = $data[4];
	$freq{$data[0]}{RA} = $data[5];
	$freq{$data[0]}{AG} = $data[6];
	$freq{$data[0]}{AA} = $data[7];
}
close ALL;

open IN,"<$need" or die($!);
while (<IN>){
	chomp;
	my @data = split /\t/,$_;
	print "$_\t$freq{$data[0]}{RG}\t$freq{$data[0]}{RA}\t$freq{$data[0]}{AG}\t$freq{$data[0]}{AA}\n";
}
close IN;
