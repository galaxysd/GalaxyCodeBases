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
my $markerdat = '/share/users/huxs/work/tiger/paper/rec.npa';
print "From [$inf] to [$inf$outf.(dat|svg)]\n";

my (@ChrOrder,%ChrIDLen,%ScaffoldLen,%ChrExists,%ChrNameLen,@ChrNameOrder);
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
	next if /^#/;
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
	$chrid = 'chr'.$chrid;
	$ChrIDLen{$items[0]} = [$chrid,$items[1]];
	$ChrNameLen{$chrid} = $items[1];
}
close I;
for (@ChrOrder) {
	print join("\t",$_,@{$ChrIDLen{$_}}),"\n";
	push @ChrNameOrder,$ChrIDLen{$_}->[0];
}
ddx \@ChrNameOrder,\%ChrNameLen;

open I,'<',$scaffnfo or die $!;
while (<I>) {
	next if /^#/;
	chomp;
	my @items = split /\t/;
	$ScaffoldLen{$items[0]} = [$items[1],$items[2]];
}
close I;
print "\n",scalar keys(%ScaffoldLen)," scaffold(s) Load.\n";

my (%Chr2Scaff,%Scaff2ChrID);
open I,'<',$inf.'.plst' or die $!;
while (<I>) {
	next if /^#/;
	chomp;
	my @items = split /\t/;
	unless (exists $ScaffoldLen{$items[3]} and exists $ChrNameLen{$items[0]}) {
		warn "[!plst]$_\n";
		next;
	}
	$Scaff2ChrID{$items[3]} = $items[0];
	my ($end,$endmin);
	if ($items[5]) {	# k
		$end = int( $items[1]-1 + $items[5]*($ScaffoldLen{$items[3]}->[0]) );
		$endmin = int( $items[1]-1 + $items[5]*($ScaffoldLen{$items[3]}->[1]) );
	} else {
		$items[5] = $items[2]?-1:1;
		$end = int( $items[1]-1 + $items[5]*($ScaffoldLen{$items[3]}->[0]) );
		$endmin = int( $items[1]-1 + $items[5]*($ScaffoldLen{$items[3]}->[1]) );
	}
	push @{$Chr2Scaff{$items[0]}},[ $items[1],$end,$endmin,$items[3],$items[5],$items[2] ];
	# start,end,endmin,scaff,k,isDiffStrand
	#push @{$RawDat{$items[0]}},[ $items[1],$end,@items[2..5],int(0.5-100*log($items[9])/log(10))/100 ];
#print join(',',@{$RawDat{$items[0]}{$items[1]}}),"\n";	# @items[2..5] -> isDiffStrand, slope(k), scaffold, pos; lg(p)
}
close I;

my %Stat;
for my $chr (keys %Chr2Scaff) {
	for my $i (0 .. $#{$Chr2Scaff{$chr}}) {
		my $arr = $Chr2Scaff{$chr}->[$i];
		++$Stat{'Total_Scaffold'};
		if ($arr->[5]) {	# isDiffStrand => $arr->[0] > $arr->[1], check with pre
			if ($i) {
				my ($pR) = sort { $b<=>$a } @{$Chr2Scaff{$chr}->[$i-1]}[0,1];
				$pR = $Chr2Scaff{$chr}->[$i-1][0];
				if ( $arr->[1] < $pR ) {
#print join(',',@{$Chr2Scaff{$chr}->[$i-1]}),",---\n";
#print join(',',@{$Chr2Scaff{$chr}->[$i]}),",$pR\n";
#$arr->[6] = '-'.$arr->[1];
#$arr->[7] = $arr->[4];
					$arr->[1] = $pR + 1;
					$arr->[4] = (1 + $arr->[1] - $arr->[0])/($ScaffoldLen{$arr->[3]}->[0]);
					if (abs($arr->[4])<0.5) {
						$arr->[1] = $arr->[2];
						$arr->[4] = (1 + $arr->[1] - $arr->[0])/($ScaffoldLen{$arr->[3]}->[0]);
						++$Stat{'Repos_Pre_Failbak'};
#$arr->[8] = 'Failbak';
					} else {
						++$Stat{'Repos_Pre'};
					}
				}
			}
		} else {	# sameStrand, check with next
			if ( $i<$#{$Chr2Scaff{$chr}} ) {
				my ($pL) = sort { $a<=>$b } @{$Chr2Scaff{$chr}->[$i+1]}[0,1];
				$pL = $Chr2Scaff{$chr}->[$i+1][0];
				if ( $arr->[1] > $pL ) {
#print join(',',@{$Chr2Scaff{$chr}->[$i+1]}),",+++\n";
#print join(',',@{$Chr2Scaff{$chr}->[$i]}),",$pL\n";
#$arr->[6] = '+'.$arr->[1];
#$arr->[7] = $arr->[4];
					$arr->[1] = $pL - 1;
					$arr->[4] = (1 + $arr->[1] - $arr->[0])/($ScaffoldLen{$arr->[3]}->[0]);
					if (abs($arr->[4])<0.5) {
						$arr->[1] = $arr->[2];
						$arr->[4] = (1 + $arr->[1] - $arr->[0])/($ScaffoldLen{$arr->[3]}->[0]);
						++$Stat{'Repos_Next_Failbak'};
#$arr->[8] = 'Failbak';
					} else {
						++$Stat{'Repos_Next'};
					}
				}
			}
		}
	}
}
ddx \%Chr2Scaff;

my %MarkerDat;
open I,'<',$markerdat or die $!;
while (<I>) {
	next if /^#/;
	chomp;
	my @items = split /\t/;
	die unless exists $ScaffoldLen{$items[1]};
	if (exists $Scaff2ChrID{$items[1]}) {
		push @{$MarkerDat{$items[1]}},[ $items[2],$items[9],int(0.5-100*log($items[9])/log(10))/100 ];	# pos,p,lg(p)
		++$Stat{'Marker_Used'};
	} else {
		++$Stat{'Marker_Skipped'};
	}
}
close I;
ddx \%MarkerDat;

# ------ BEGIN PLOT --------
my $Xrange = 2000;
my $Yrange = 100;
my $ArrowLen = 20;
my $OutBorder = 24;
my $InBorder = 36;

my %PlotDat;
for my $chr (keys %Chr2Scaff) {
	for my $items (@{$Chr2Scaff{$chr}}) {	# start,end,endmin,scaff,k,isDiffStrand
		die unless exists $MarkerDat{$$items[3]};
		my $x;
		push $PlotDat{$chr}{$$items[3]},[];
	}
}
ddx \%Stat;

__END__
