#!/usr/bin/env perl
use strict;

my $input = shift;
my $output = shift;

open IN,"sort -Vk2 $input|" or die($!);
open OUT,">$output" or die($!);
print OUT "#CHROM, START0, END, Name, Score\n";
while (<IN>){
	chomp;
	my @data = split /\t/,$_;
	if ($#data == 3 && $data[$#data] > 0){
		print OUT "$data[1]\t",$data[2] - 1,"\t$data[3]\t$data[0]\n";
	}else{
		print OUT "$data[1]\t",$data[2] - 1,"\t$data[2]\t$data[0]\n";
	}
}
close IN;
close OUT;

system("bgzip $output && tabix $output.gz");
