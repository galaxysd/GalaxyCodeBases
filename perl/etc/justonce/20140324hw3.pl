#!/usr/bin/perl -w
use strict;
use warnings;
use Data::Dump;

die "Usage: $0 <s> <t>\nW(AA), W(Aa) and W(aa) is 1-s, 1, 1-t respectively.\n" if @ARGV<2;

my ($s,$t)=@ARGV;

my ($w11,$w12,$w22) = (1-$s, 1, 1-$t);
my ($AA,$Aa,$aa) = (0.50, 0.40, 0.10);

my $p = $AA + $Aa/2;
my $q = $aa + $Aa/2;

my $avgW = $p*$p*$w11 + 2*$p*$q*$w12 + $q*$q*$w22;

print "GT: ($AA,$Aa,$aa) => ($p,$q)\nFit: ($s,$t) => ($w11,$w12,$w22), $avgW\n";

for my $i ( 1 .. 2 ) {
	my $avgW = $p*$p*$w11 + 2*$p*$q*$w12 + $q*$q*$w22;
	my $np = ( $p*$p*$w11 + $p*$q*$w12 ) / $avgW;
	my $nq = ( $p*$q*$w12 + $q*$q*$w22 ) / $avgW;
	$AA = $p*$p; $Aa = 2*$p*$q; $aa = $q*$q;
	print "$i\t-> ($np, $nq) ",$np+$nq,"\t($AA,$Aa,$aa)\tavg(w)=$avgW\n";
	$p = $np; $q = $nq;
}
