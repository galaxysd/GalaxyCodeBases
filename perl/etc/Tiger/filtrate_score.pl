#!/usr/bin/perl
use strict;
use warnings;

while (<>) {
	chomp;
	my @a = split /\t/;
	print "$_\n" if $a[-1] > 3000;
}
