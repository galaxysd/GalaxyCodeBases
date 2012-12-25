#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Galaxy::IO;

die "Usage: $0 <pad to len> <fq.gz> [out midfix]\n" if @ARGV<2;
my $padto=shift;
my $inf=shift;
my $outf=shift;

my $main = $inf;
$main =~ s/\.fq(\.gz)$//;
if ($outf) {
	$main .= ".$outf"
}
$main .= ".$padto.fq.gz";
warn "\@$padto. From [$inf] to [$main]\n";

my $infh = openfile($inf);
open OUT,'|-',"gzip -9c > $main" or die "Error opening $main: $!\n";
my ($a,$b,$c,$d);
while (<$infh>) {
	$a=$_;
	chomp($b=<$infh>);
	$c=<$infh>;
	chomp($d=<$infh>);
	if (length($b) < $padto) {
		my $t = $padto - length($b);
		$b .= 'N' x $t;
		$d .= '#' x $t;
	} elsif (length($b) == $padto) {
		;
	} else {die;}
	print OUT "$a$b\n$c$d\n";
}

close $infh;
close OUT;
