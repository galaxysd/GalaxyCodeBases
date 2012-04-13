#!/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120412
=cut
use strict;
use warnings;
#use Data::Dump qw(ddx);

die "Usage: $0 <Sperm count> <Site count>\n" if @ARGV<2;
my $indi=shift;
my $site=shift;

my $RateCMperMb=2.1;
my $ChrLen=249250621;

my $Rate = ($ChrLen/$site/1000000)*$RateCMperMb/100;
my (@Indi,%S);
open O,'>',"gt_${indi}_${site}.dat" or die;
print O "#Rate=$Rate, Distence=",int($ChrLen/$site),"\nRecombined\tSNP\t\@Genotype_of_${indi}\n";
push @{$Indi[0]},0 for (0 .. $indi-1);
for my $snp (1..$site) {
	my @t=();
	my $sumC=0;
	for my $i (0 .. $indi-1) {
		my $x=(rand(1)>$Rate)?0:1;
		if ($x) {	# 1:recombine
			push @t,($Indi[$snp-1][$i])?0:1;
		} else {
			push @t,($Indi[$snp-1][$i])?1:0;
		}
		$sumC += $x;
	}
	print O join("\t",$sumC,$snp,@t),"\n";
	++$S{$sumC};
	$Indi[$snp]=\@t;
}
close O;

open O,'>',"gt_${indi}_${site}.stat" or die;
print O "Recombine Times:\n";
for my $c (sort {$a<=>$b} keys %S) {
	print O "$c\t$S{$c}\n";
}

my %Sgt;
print O "\nGenoTypes:\n";
for my $snp (1..$site) {
	for my $i (0 .. $indi-1) {
		++$Sgt{$Indi[$snp][$i]};
	}
}
print O "0: $Sgt{0}, 1: $Sgt{1}. Ratio: ",$Sgt{0}/$Sgt{1},"\n";
close O;

