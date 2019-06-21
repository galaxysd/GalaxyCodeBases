#!/usr/bin/env perl
use strict;
use warnings;

use Data::Dump qw(ddx);

die "Usage: $0 <in.tsv>\n" if @ARGV<1;
#my ($inf) = @ARGV;
while (<>) {
	next if /^#/;
	chomp;
	my @d = split /\t/;
	my $freq2 = pop @d;
	my $Alts = pop @d;
	my @AltGTs = split /,/,$Alts;
	my %GTfreq;
	$GTfreq{$AltGTs[0]} = $freq2;
	my $sum2 = $freq2;
	if (@AltGTs > 1) {
		for (1 .. $#AltGTs) {
			$GTfreq{$AltGTs[$_]} = 0.01;
			$sum2 += 0.01;
		}
	}
	$GTfreq{$d[3]} = 1 - $sum2;
	my @allels = sort { $GTfreq{$a} <=> $GTfreq{$b} } keys %GTfreq;
	my @af = map {"$_\t$GTfreq{$_}"} @allels;
	$d[0] = 'chr' . $d[0];
	print join("\t",@d[2,0,1],@af),"\n";
}

__END__
./cpe.pl freq/MGIEasy.freq.txt >freq/MGIEasy.freq.tsv
