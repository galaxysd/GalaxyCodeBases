#!/usr/bin/perl
use strict;
use warnings;
use Text::CSV;
use Data::Dump qw(ddx);

my $glist = 'input.hg38_multianno.txt';
my $finput = 'ores.csv';
my $TheDis = 10000;

my $csv = Text::CSV->new ({ sep => "\t",binary => 1 }) or die "Cannot use CSV: ".Text::CSV->error_diag ();
my %iDat;

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
		#my @dat = ($row->[0],$row->[1],[@to]);
		#ddx \@dat;
		$iDat{$row->[0]}{$row->[1]} = [@to];
	}
}
close $fh;
#ddx \%iDat;

$csv->sep (",");

open OI,'>',"res.witnin" or die $!;
open OO,'>',"res.outside" or die $!;

open $fh, "<:encoding(utf8)",$finput or die "$!\n";
@cols = @{$csv->getline ($fh)};
#$csv->column_names (@cols);
while ( my $row = $csv->getline($fh) ) {
	#ddx $row;
	if (exists $iDat{$row->[1]}) {
		my $flag = 0;
		my $p1 = $row->[3] - $TheDis;
		my $p2 = $row->[4] + $TheDis;
		$p1 = 0 if $p1 < 0;
		for my $p ($p1 .. $p2) {	# I am in a hurry.
			if (exists $iDat{$row->[1]}{$p}) {
				my $t = join(",",@{$iDat{$row->[1]}{$p}});
				print OI join("\t",$t,$p,@$row),"\n";
				++$flag;
			}
		}
		unless ($flag) {
			print OO join("\t",'NA',0,@$row),"\n";
		}
	}
}
close $fh;
close OI;
close OO;
