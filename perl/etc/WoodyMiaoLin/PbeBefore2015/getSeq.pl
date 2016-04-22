#!/usr/bin/perl
use strict;
use warnings;

open I, "<", $ARGV[0];
#open O1, ">", $ARGV[1];
#open O2, ">", $ARGV[2];
open O3, ">", $ARGV[1];

my %seq;
while (<I>) {
	chomp;
	s/>//;
	$seq{$_} = <I>;
}
close I;

#my @north = qw/ Mir4 Pbe106 Pbe107 Pbe108 Pbe11 Pbe111 Pbe114 Pbe115 Pbe116 Pbe117 Pbe118 Pbe119 Pbe120 Pbe123 Pbe124 Pbe125 Pbe126 Pbe127 Pbe128 Pbe129 Pbe13 Pbe130 Pbe131 Pbe132 Pbe133 Pbe134 Pbe135 Pbe136 Pbe137 Pbe138 Pbe139 Pbe14 Pbe140 Pbe142 Pbe144 Pbe145 Pbe146 Pbe148 Pbe149 Pbe15 Pbe152 Pbe153 Pbe154 Pbe155 Pbe156 Pbe157 Pbe158 Pbe159 Pbe16 Pbe160 Pbe161 Pbe166 Pbe167 Pbe169 Pbe170 Pbe171 Pbe172 Pbe173 Pbe174 Pbe175 Pbe176 Pbe179 Pbe180 Pbe183 Pbe187 Pbe188 Pbe191 Pbe193 Pbe2 Pbe24 Pbe25 Pbe26 Pbe27 Pbe28 Pbe29 Pbe3 Pbe30 Pbe31 Pbe32 Pbe34 Pbe35 Pbe37 Pbe38 Pbe39 Pbe40 Pbe41 Pbe42 Pbe43 Pbe51 Pbe52 Pbe53 Pbe54 Pbe55 Pbe56 Pbe57 Pbe58 Pbe60 Pbe61 Pbe62 Pbe63 Pbe64 Pbe65 Pbe66 Pbe67 Pbe68 Pbe69 Pbe7 Pbe70 Pbe71 Pbe72 Pbe73 Pbe74 Pbe75 Pbe76 Pbe78 Pbe79 Pbe8 Pbe82 Pbe89 Pbe91 Pbe92 Pbe93 Pbe94 Pbe95 /;
#my @sunda = qw/ Pbe100 Pbe101 Pbe102 Pbe113 Pbe17 Pbe18 Pbe198 Pbe20 Pbe204 Pbe205 Pbe21 Pbe44 Pbe48 Pbe50 Pbe84 Pbe85 Pbe86 Pbe88 Pbe97 Pbe99 /;
#my @amur = qw/ Mir4 Pbe123 Pbe157 Pbe158 Pbe159 Pbe160 Pbe161 Pbe175 Pbe24 Pbe25 Pbe26 Pbe27 Pbe28 Pbe29 Pbe30 Pbe31 Pbe32 Pbe33 Pbe34 Pbe35 Pbe40 Pbe42 Pbe43 Pbe54 Pbe55 Pbe91 Pbe92 Pbe93 Pbe94 Pbe95 /;
my @indochina = qw/ Pbe9 Pbe11 Pbe13 Pbe130 Pbe139 Pbe14 Pbe15 Pbe37 Pbe39 Pbe51 Pbe56 Pbe7 Pbe76 Pbe8 Pbe89 Pbe60 Pbe61 Pbe71 Pbe119 Pbe155 Pbe68 Pbe69 Pbe153 Pbe156 Pbe166 Pbe171 Pbe10 Pbe6 Pbe16 Pbe52 Pbe82 Pbe135 Pbe2 Pbe3 Pbe70 Pbe74 Pbe174 Pbe114 Pbe116 Pbe67 Pbe132 Pbe62 Pbe63 Pbe64 Pbe167 Pbe169 Pbe170 Pbe172 Pbe173 Pbe115 Pbe118 Pbe120 Pbe152 Pbe72 Pbe78 Pbe79 Pbe131 Pbe133 Pbe112 Pbe108 Pbe111 Pbe106 Pbe107 Pbe140 Pbe146 Pbe117 Pbe142 Pbe38 Pbe134 Pbe136 Pbe145 Pbe66 Pbe127 Pbe129 Pbe125 Pbe126 Pbe128 Pbe124 Pbe41 Pbe65 Pbe137 Pbe138 Pbe144 Pbe148 Pbe58 Pbe141 Pbe149 Pbe73 Pbe75 Pbe90 Pbe59 /;

#foreach (sort @north) {
#	print O1 ">$_\n$seq{$_}" if $seq{$_};
#}
#foreach (sort @sunda) {
#	print O2 ">$_\n$seq{$_}" if $seq{$_};
#}

foreach (sort @indochina) {
	print O3 ">$_\n$seq{$_}" if $seq{$_};
}

#close O1;
#close O2;
close O3;
