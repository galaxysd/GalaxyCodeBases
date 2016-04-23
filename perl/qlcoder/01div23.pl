#!/usr/bin/env perl
use strict;
use warnings;

my ($i,$z)=(0,0);

while ($i <= 2333) {
	++$z;
	if (( int($z/2) == $z/2 ) or ( int($z/3) == $z/3 )) {
		++$i;
		print "$i: $z\n";
	}
}
