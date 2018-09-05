#!/usr/bin/perl
use strict;
use warnings;
use Text::CSV;
use Data::Dump qw(ddx);

my $glist = 'input.hg38_multianno.txt';
my $finput = 'ores.csv';
my $TheDis = 10000;

my $csv = Text::CSV->new ({ sep => "\t",binary => 1 }) or die "Cannot use CSV: ".Text::CSV->error_diag ();

open my $fh, "<:encoding(utf8)", $glist or die "$!\n";
my @cols = @{$csv->getline ($fh)};
while ( my $row = $csv->getline($fh) ) {
	my @t = split /,/,$row->[6];
	my @to;
	for my $i (@t) {
		next if $i eq 'NONE(dist=NONE)';
		$i =~ /^([\w\-\.]+)(\(dist=(\d+)\))?$/ or die $i;
		if ($2?($3 <= $TheDis):1) {
			push @to,$1;
		}
	}
	if (scalar @to > 0) {
		my @dat = ($row->[0],$row->[1],[@to]);
		ddx \@dat;
	}
}
close $fh;

