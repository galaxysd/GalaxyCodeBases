#!/usr/bin/env perl
use strict;
use warnings;

use IPC::Open2;
my  $pid = open2( \*READER, \*WRITER, "./alnmethly" );
    WRITER->autoflush(); # default here, actually
my $got = "";
my $input = " ";

print WRITER "a,gtttaacccgg\n,gttCggttCg\nxxx\n";

while(<READER>) {
	chomp;
	if (/^Path(\d): ([IDMmR]+),(\d+)$/) {
		print "[$_] [$1] [$2] [$3]\n";
	}
}
