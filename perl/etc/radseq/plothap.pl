#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Data::Dump qw(ddx);
use Galaxy;

die "Usage: $0 <marker dat> <tassel dump> [out midfix]\n" if @ARGV<2;
my $markf=shift;
my $inf=shift;
my $outf=shift;

if ($outf) {
	$outf = ".$outf";
} else {
	$outf = '';
}

#my $scaffnfo = '/bak/seqdata/2012/tiger/120512_TigerRefGenome/chr.nfo';
my $scaf2locuslst = '/share/users/huxs/work/tiger/20120910/tighap.lst';
print "From [$markf][$inf] to [$inf$outf.(dat|svg)]\n";

my %Stat;

my %LD;
open I,'<',$inf or die $!;
while (<I>) {
	next if /^Locus1\t/;
	chomp;
	my @items = split /\t/;
	die scalar @items if @items != 17;
	my ($chr1,$pos1,$id1,undef,undef,undef,$chr2,$pos2,$id2,undef,undef,undef,$dist,$rsq,$dp,$p,$N) = @items;
	my $x = [$chr1,$pos1];
	my $y = [$chr2,$pos2];
#ddx [$x,$y];
	($x,$y) = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } ($x,$y);
#ddx [$x,$y];
	$LD{$$x[0]}{$$x[1]}{$$y[0]}{$$y[1]} = [$dist,$rsq,$dp,$p,$N];
}
ddx \%LD;

__END__



