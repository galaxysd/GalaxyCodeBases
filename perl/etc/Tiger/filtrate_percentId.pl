#!/usr/bin/perl
use strict;
use warnings;

while (<>) {
	my @a = split /\t/;
	print $_ if $a[2] > 97;
}
