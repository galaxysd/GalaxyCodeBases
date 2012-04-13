#!/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120412
=cut
use strict;
use warnings;
#use Data::Dump qw(ddx);

die "Usage: $0 <in file> <out file>\n" if @ARGV<2;
my $in=shift;
my $out=shift;

my @Indi;
my ($indi,$site)=(0,0);

open I,'<',$in or die;
while (<I>) {
	next if /^#/;
	chomp;
	my (undef,undef,$id,@GT)=split /\t/;
	$Indi[$id] = \@GT;
	die if $indi && $indi != @GT;
	$indi = scalar @GT;
}
$site = $#Indi;
push @{$Indi[0]},0 for (0 .. $indi-1);
print "Loaded: $site SNP sites * $indi Sperms.\n";

