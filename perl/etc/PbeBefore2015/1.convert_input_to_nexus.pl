#!/usr/bin/perl
use strict;
use warnings;

open I, "<", "../f.Gronau_2011_analysis/4.G-PhoCS/sequence.input";
open O, ">", "sequence.nex";

print O "#NEXUS
BEGIN DATA;
Dimensions ntax=4 nchar=32648000;
Format datatype=DNA interleave=yes;
Matrix\n";

while (<I>) {
	next unless /^[A-Z]|^\n/;
	print O $_;
}

print O ";\nEND;\n";

close I;
close O;
