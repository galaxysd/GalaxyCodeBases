#!/bin/env perl
use strict;
use warnings;
use IO::Unread qw(unread);
use Data::Dump qw(ddx);

#die "Usage: $0 <input1> <input2> <output>\n" if @ARGV < 3;
#my ($inf1,$inf2,$outf)=@ARGV;

my (%Files,%Fname2ID);
open L,'<','kmerfreq.lst' or die $!;
while (<L>) {
	chomp;
	my @t = split /\//,$_;
	my @k = split /\./,$t[2];
	if ($t[2] =~ /^Donor\./) {
		push @{$Files{$k[1]}{'Donor'}},$_;
	} else {
		push @{$Files{$k[1]}{$t[1]}},$_;
	}
	$Fname2ID{$_} = $k[0];
}
ddx \%Files;
close L;

my $prefix = 'perl plotkmercnt.pl 1e-8 ';
my $outpath = 'out';
my $linending = " &\n";

sub getout($$$$) {
	my ($x,$y,$k,$outpath) = @_;
	return join(' ',$x,$y,"$outpath/$k-$Fname2ID{$x}-$Fname2ID{$y}");
}

open OUT,'>','run.sh' or die $!;
for my $k ( keys %Files ) {
	my @types = keys %{$Files{$k}};
	for my $i ( 0 .. ($#types-1) ) {
		for my $j ( ($i+1) .. $#types ) {
			print "$k: $i,$j\n";
			my $arrA = $Files{$k}->{$types[$i]};
			my $arrB = $Files{$k}->{$types[$j]};
			# pairwise
			for my $x ( 0 .. $#$arrA ) {
				for my $y ( 0 .. $#$arrB ) {
					print "p $arrA->[$x], $arrB->[$y]\n";
					print OUT $prefix,getout($arrA->[$x],$arrB->[$y],$k,$outpath),$linending;
				}
			}
			# within A
			for my $x ( 0 .. ($#$arrA-1) ) {
				for my $y ( ($x+1) .. $#$arrA ) {
					print "wA $arrA->[$x], $arrA->[$y]\n";
					print OUT $prefix,getout($arrA->[$x],$arrA->[$y],$k,$outpath),$linending;
				}
			}
			# within B
			for my $x ( 0 .. ($#$arrB-1) ) {
				for my $y ( ($x+1) .. $#$arrB ) {
					print "wB $arrB->[$x], $arrB->[$y]\n";
					print OUT $prefix,getout($arrB->[$x],$arrB->[$y],$k,$outpath),$linending;
				}
			}
		}
	}
}
close OUT;
