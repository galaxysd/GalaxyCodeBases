#!/usr/bin/perl
use strict;
use warnings;

die("\nUsage: $0 <input_f> <output.svg>\n\n") if (@ARGV < 2);

my $in = shift;
my $out = shift;
open I, "<", $in;
while (<I>) {
	last if /\(%Miss\) :  Inferred clusters/;
}

my @color = ("#FF0000", "#00FF00", "#0000FF", "#00FFFF", "#FF00FF", "#FFFF00", "#800000", "#008000", "#000080", "#008080", "#800080", "#808000");
my @PO; # Probability of inferred clusters with original order
my @pg; # Probability of inferred clusters grouped by numerical order of the max probability
my $num_clust; # Number of clusters
my $num_indiv; # Number of individuals
while (<I>) {
	chomp;
	last unless $_;
	my @a = split /:  /;
	$a[1] =~ s/ $//;
	my @b = split / /, $a[1];
	$num_clust = @b unless $num_clust;
	++$num_indiv;
	unshift @b, $num_indiv;
	my $max_order = 1; # The numerical order of the max probability
	foreach (1 .. $num_clust) {
		if ($b[$_] > $b[$max_order]) {
			$max_order = $_;
		}
	}
	push @{$pg[$max_order]}, \@b;
	push @PO, \@b;
}
close I;

my @PQ; # Probability of inferred clusters sort by Q
foreach (1 .. $num_clust) {
	next unless $pg[$_];
	push @PQ, sort {${$b}[$_] <=> ${$a}[$_]} @{$pg[$_]};
}

$in =~ s/^.*\///;
open O, ">", $out;
# SVG attribute definitions
my $w = 10 * $num_indiv + 20;
my $h = 240;
print O "<?xml version=\"1.0\"?>\n";
print O "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"$w\" height=\"$h\">\n";
print O "<g transform=\"translate(10,10)\">\n\n";

# Print histogram
foreach my $indiv (0 .. $num_indiv-1) {
	my $xr = 10*($indiv);
	my $yrPO = 0;
	my $yrPQ = 115;
	my $xt = $xr + 1;
	foreach my $clust (1 .. $num_clust) {
		my $hPO = 100*$PO[$indiv][$clust];
		my $hPQ = 100*$PQ[$indiv][$clust];
		print O "<rect x=\"$xr\" y=\"$yrPO\" width=\"10\" height=\"$hPO\" fill=\"$color[$clust-1]\"/>\n";
		print O "<rect x=\"$xr\" y=\"$yrPQ\" width=\"10\" height=\"$hPQ\" fill=\"$color[$clust-1]\"/>\n";
		$yrPO += $hPO;
		$yrPQ += $hPQ
	}
	my $numPO = sprintf("%03d", $PO[$indiv][0]);
	my $numPQ = sprintf("%03d", $PQ[$indiv][0]);
	print O "<text x=\"$xt\" y=\"105\" font-family=\"Courier\" font-size=\"4\" fill=\"black\">$numPO</text>\n";
	print O "<text x=\"$xt\" y=\"220\" font-family=\"Courier\" font-size=\"4\" fill=\"black\">$numPQ</text>\n";
}

print O "\n</g>\n</svg>\n";
close O;
