#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Data::Dump qw(ddx);

die "Usage: $0 <prefix> [out midfix]\n" if @ARGV<1;
my $inf=shift;
my $outf=shift;

if ($outf) {
	$outf = ".$outf";
} else {
	$outf = '';
}

my $scaffnfo = '/bak/seqdata/2012/tiger/120512_TigerRefGenome/chr.nfo';
print "From [$inf] to [$inf$outf.(dat|svg)]\n";

my (@ChrOrder,%ChrIDLen,%SacffoldLen,%ChrExists);
open I,'<',$inf.'.chrorder' or die $!;
while (<I>) {
	next if /^#/;
	chomp;
	push @ChrOrder,$_;
	++$ChrExists{$_};
}
close I;
open I,'<',$inf.'.chrnfo' or die $!;
while (<I>) {
	chomp;
	my @items = split /\t/;
	next unless exists $ChrExists{$items[0]};
	my $chrid;
	if ($items[9] =~ /chromosome (\w+)/) {
		$chrid = $1;
	} elsif ($items[9] =~ /mitochondrion/) {
		$chrid = 'M';
	} elsif ($items[0] =~ /^gi\|\d+\|\w+\|([^|]+)\|/) {
		$chrid = $1;
	} else {
		$chrid = $items[0];
	}
	$ChrIDLen{$items[0]} = [$items[1],$chrid];
}
close I;
print join("\t",$_,@{$ChrIDLen{$_}}),"\n" for (@ChrOrder);

