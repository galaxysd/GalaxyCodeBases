#!/usr/bin/perl
use strict;
use warnings;

open IP, "<", "phase.ht";
open IM, "<", "manual.ht";

my (@id, @phase, @manual, @compare, @conflict, @out);

my $n = 0;
while (<IP>) {
	chomp;
	s/^0 //;
	push @id, $_;
	my $h1 = <IP>;
	my $h2 = <IP>;
	chomp ($h1, $h2);
	$phase[$n][0] = [split / /, $h1];
	$phase[$n][1] = [split / /, $h2];
	++$n;
}
close IP;

while (<IM>) {
	chomp;
	my @a = split /\t/;
	my $n = 0;
	foreach (@a) {
		if ($n == int($n)) {
			push @{$manual[int($n)][0]}, $_;
		} else {
			push @{$manual[int($n)][1]}, $_;
		}
		$n += 0.5;
	}
}
close IM;

foreach my $n (0 .. 11) {
	$conflict[$n][0] = 0;
	$conflict[$n][1] = 0;
	foreach my $p (0 .. 130) {
		foreach my $a1 (0 .. 1) {
			foreach my $a2 (0 .. 1) {
				if ($phase[$n][$a1][$p] eq $manual[$n][$a2][$p]) {
					$compare[$n]{"c$a1$a2"}[$p] = 0;
				} else {
					$compare[$n]{"c$a1$a2"}[$p] = 1;
				}
			}
		}
		$conflict[$n][0] += $compare[$n]{c00}[$p] + $compare[$n]{c11}[$p];
		$conflict[$n][1] += $compare[$n]{c01}[$p] + $compare[$n]{c10}[$p];
	}
	if ($conflict[$n][1] > $conflict[$n][0]) {
		$out[$n][0][0] = $phase[$n][0];
		$out[$n][0][1] = $manual[$n][0];
		$out[$n][0][2] = $compare[$n]{c00};
		$out[$n][1][0] = $phase[$n][1];
		$out[$n][1][1] = $manual[$n][1];
		$out[$n][1][2] = $compare[$n]{c11};
	} else {
		$out[$n][0][0] = $phase[$n][0];
		$out[$n][0][1] = $manual[$n][1];
		$out[$n][0][2] = $compare[$n]{c01};
		$out[$n][1][0] = $phase[$n][1];
		$out[$n][1][1] = $manual[$n][0];
		$out[$n][1][2] = $compare[$n]{c10};
	}
}

open O, ">", "compare.ht";

foreach (0 .. 11) {
	print O "$id[$_]\n";
	print O "hyplotype0\tphase\t", (join " ", @{$out[$_][0][0]}), "\n";
	print O "hyplotype0\tmanual\t", (join " ", @{$out[$_][0][1]}), "\n";
	print O "hyplotype0\tconfl.\t", (join " ", @{$out[$_][0][2]}), "\n\n";
	print O "hyplotype1\tphase\t", (join " ", @{$out[$_][1][0]}), "\n";
	print O "hyplotype1\tmanual\t", (join " ", @{$out[$_][1][1]}), "\n";
	print O "hyplotype1\tconfl.\t", (join " ", @{$out[$_][1][2]}), "\n\n";
}
close O;
