#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

die "Usage: $0 <depth file>\n" if @ARGV <1;

my $in = shift;

open I,'<',$in or die;
my (%Depdat,%Ranges);
while (<I>) {
	chomp;
	my ($chr,$pos,$depth) = split /\t/;
	$Depdat{$chr}{$pos} = $depth;
}
close I;
for my $chr (sort keys %Depdat) {
	my @poses = sort {$a <=> $b} keys %{$Depdat{$chr}};
	my ($st,$ed);
	while (@poses) {
		$st = shift @poses;
		$ed = shift @poses;
		while (@poses && $poses[0] - $ed == 1) {
			$ed = shift @poses;
		}
		push @{$Ranges{$chr}},[$st,$ed];
	}
}
ddx \%Ranges;
