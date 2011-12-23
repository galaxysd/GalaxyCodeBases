#!/bin/env perl
use strict;
use warnings;

die "Usage: $0 <time ??m?.?s>\n" if @ARGV < 1;

my $totalsec=0;

while($_=shift @ARGV) {
	print "[$_]: ";
	/((\d+)m)?([.\d]+)s/ or (print "NULL\n" and next);
	#print "$1, $2, $3\n";
	my $sec;
	if (defined $1) {
		$sec += $2 * 60;
	}
	$sec += $3;
	print "$sec s\n";
	$totalsec += $sec;
}
print "Total: $totalsec s\n";

