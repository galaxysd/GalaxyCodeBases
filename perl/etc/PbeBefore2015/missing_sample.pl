#!/usr/bin/perl
use strict;
use warnings;

open I, "<", "CALB-CHRNA-CMA-GNB.fasta";
open O1, ">", "CALB-CHRNA-CMA-GNB_missMore100.fasta";
open O2, ">", "CALB-CHRNA-CMA-GNB_missLess100.fasta";

while (<I>) {
	my $seq = <I>;
	my $miss = @{[$seq =~ /-/g]};
	if ($miss > 100) {
		print O1 $_, $seq;
	} else {
		print O2 $_, $seq;
	}
}
close I;
close O1;
close O2;
