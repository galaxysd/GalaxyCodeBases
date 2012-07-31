#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;

die "Usage: $0 <ec list> <out>\n" if @ARGV<1;
my $eclst=shift;
my $outfs=shift;

open O,'>',$outfs or die "Error opening $outfs : $!\n";

my (%Markers,@ChrOrder);
open L,'<',$eclst or die;
while (<L>) {
	if (/^(#|$)/) {
		print O $_;
		next;
	}
	my ($chr,$pos,$strand,$mark,$count,$sampleCNT) = split /\t/;
	push @ChrOrder,$chr unless exists $Markers{$chr};
	push @{$Markers{$chr}{$pos}},[$strand,$mark,$count,$sampleCNT];
}
close L;

for my $chr (@ChrOrder) {
	my ($Fstrand,$Fmark,$Fcount,$FsampleCNT);
	for my $pos (keys %{$Markers{$chr}}) {
		next if @{$Markers{$chr}{$pos}} == 1;
		($Fstrand,$Fmark,$Fcount,$FsampleCNT)=(0,0,0,0);
		for (@{$Markers{$chr}{$pos}}) {
			my ($strand,$mark,$count,$sampleCNT)=@$_;
			if ($strand eq '+') {
				$Fstrand |= 1;
			} elsif ($strand eq '-') {
				$Fstrand |= 2;
			} else {
				$Fstrand = 3;
				$mark -= 10000;
			}
			$Fmark += $mark;
			$Fcount += $count;
			$FsampleCNT = $sampleCNT if $FsampleCNT < $sampleCNT;
		}
		if ($Fstrand == 1) {
			$Fstrand = '+';
		} elsif ($Fstrand == 2) {
			$Fstrand = '-';
		} else {
			$Fstrand = 'A';
			$Fmark += 10000;
		}
		$Markers{$chr}{$pos} = [[$Fstrand,$Fmark,$Fcount,$FsampleCNT]];
	}
}

for my $chr (@ChrOrder) {
	for my $pos (sort {$a <=> $b} keys %{$Markers{$chr}}) {
		print O join("\t",$chr,$pos,@{$Markers{$chr}{$pos}->[0]}),"\n";
	}
}
close O;
