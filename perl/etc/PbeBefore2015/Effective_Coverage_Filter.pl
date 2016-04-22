#!/usr/bin/perl
use strict;
use warnings;

my $in = shift;
my $out = shift;

open I, "-|", "gzip -dc $in";
open O, "|-", "gzip -9c >$out";

my $head = <I>;
print O $head;
while (<I>) {
	chomp;
	my @a = split /\t/;
	next if $a[-4] < 5; # Next if NumGoodReads < 5
	my @rq = split //, $a[-2]; # Read Quality Scores
	my @aq = split //, $a[-1]; # Align Quality Scores
	my $nr = @rq; # Number of Reads
	my $ec; # Effective Coverage
	foreach (0 .. $nr-1) {
		my $rq = ord($rq[$_]) - 33;
		my $aq = ord($aq[$_]) - 33;
		my $c = (1 - (10**(-$rq/10))) * (1 - (10**(-$aq/10)));
		$ec += $c;
	}
	if ($ec < 5) {
		next;
	} else {
		print O "$_\n";
	}
}
close I;
close O;
