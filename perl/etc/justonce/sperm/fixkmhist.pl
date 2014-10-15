#!/bin/env perl
use strict;
use warnings;

die "Usage: $0 <input>\n" if @ARGV < 1;
my ($inf)=@ARGV;

open IN,'<',$inf or die "Error opening $inf: $!\n";
open OUTF,'>',$inf.'.fixed' or die "Error opening $inf.fixed: $!\n";

my ($HistSum,$Check1,$Check2)=(0,0,0);
while (<IN>) {
	if (/^#/) {
		print OUTF $_;
		if (/\bCntOfFoundKmer\t(\d+)\t/) {
			$HistSum = $1;
		}
	} else {
		chomp;
		my @DAT = split /\t/;
		push @DAT,$DAT[2];
		$DAT[2] = $DAT[1]/$HistSum;
		$Check1 += $DAT[2];
		$Check2 += $DAT[4];
		$DAT[3] = $Check1;
		print OUTF join("\t",@DAT),"\n";
	}
}
die if $HistSum == 0;
print "$HistSum,$Check1,$Check2\n";
close IN;
close OUTF;

__END__
find *.hist | xargs -n1 ./fixkmhist.pl
