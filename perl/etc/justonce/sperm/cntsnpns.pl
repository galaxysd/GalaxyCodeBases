#!/usr/bin/env perl
use strict;
use warnings;

die "Usage: $0 <stat.txt> <flag>\n" if @ARGV < 2;
my ($inf,$flag)=@ARGV;

my $sum=0;

open I,'<',$inf or die $!;
while (<I>) {
	next if /\[/;
	chomp;
	my @a = split /\t/;
	if ( ($a[2] & $flag)==$flag ) {
		$sum += $a[1];
		print join("\t",$a[0],$a[2],$a[1],$sum),"\n";
	}
}
close I;
print "Sum = $sum\n";
