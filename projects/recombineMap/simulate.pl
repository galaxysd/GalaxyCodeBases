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
my (@Indi,%S,@Parent);
open O,'>',"gt_${indi}_${site}.gt" or die;
open OR,'>',"gt_${indi}_${site}.dat" or die;
print O "#Rate=$Rate, Distence=",int($ChrLen/$site),"\n#recombRate\tRecombined\tSNPid\t\@Genotype_of_${indi}\n";
print OR "#Rate=$Rate, Distence=",int($ChrLen/$site),"\n#Parent\tRecombined\tSNPid\t\@Genotype_of_${indi}\n";
push @{$Indi[0]},0 for (0 .. $indi-1);
my $NewRate = $Rate;
my $SumRate=0;
for my $snp (1..$site) {
	my @t=();
	my $sumC=0;
	$Parent[$snp]=int(rand(2));
	$NewRate = ($Rate*$snp)-$SumRate;	# a feedback is better, EP(?)
	for my $i (0 .. $indi-1) {
		my $x=(rand(1)>$NewRate)?0:1;
		if ($x) {	# 1:recombine
			push @t,($Indi[$snp-1][$i])?0:1;
		} else {
			push @t,($Indi[$snp-1][$i])?1:0;
		}
		$sumC += $x;
	}
	$SumRate += $sumC/$indi;
	#$NewRate = ($Rate*($snp+1))-$SumRat;
	print O join("\t",sprintf("% .6f",$NewRate),$sumC,$snp,@t),"\n";
	++$S{$sumC};
	$Indi[$snp]=\@t;
	print OR join("\t",$Parent[$snp],$sumC,$snp,map {$_^$Parent[$snp]} @t),"\n";
}
close O;
close OR;

my $t=0;
my $cnt=$indi*$site;
open O,'>',"gt_${indi}_${site}.stat" or die;
print O "Rate=$Rate, Distence=",int($ChrLen/$site),"\ncM pre Mbp:",$RateCMperMb,"\nRecombine Times and Rate:\n";
for my $c (sort {$a<=>$b} keys %S) {
	print O "$c\t$S{$c}\t",$c/$indi,"\n";
	$t += $c*$S{$c};
}
print O "Average Rate: ",$t/$cnt,"\n";
print O "Average Recombine Rate: ",100*$t/$indi/($ChrLen/1000000)," (cM pre Mbp)\n\n";

my %Sgt;
print O "GenoTypes:\n";
for my $snp (1..$site) {
	for my $i (0 .. $indi-1) {
		++$Sgt{$Indi[$snp][$i]};
	}
}
print O "0: $Sgt{0}, 1: $Sgt{1}.\nRatio: ",$Sgt{0}/$Sgt{1},"\n";
close O;
