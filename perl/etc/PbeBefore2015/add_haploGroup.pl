#!/usr/bin/perl
use strict;
use warnings;

open I, "<", "a.txt";
open O, ">", "a1.txt";

my %group = (
"PbeMT1"	,	"N1"	,
"PbeMT2"	,	"N1"	,
"PbeMT3"	,	"N1"	,
"PbeMT4"	,	"N1"	,
"PbeMT5"	,	"N1"	,
"PbeMT6"	,	"N1"	,
"PbeMT7"	,	"N1"	,
"PbeMT8"	,	"N1"	,
"PbeMT9"	,	"N1"	,
"PbeMT10"	,	"N1"	,
"PbeMT11"	,	"N1"	,
"PbeMT12"	,	"N1"	,
"PbeMT13"	,	"N1"	,
"PbeMT14"	,	"N1"	,
"PbeMT15"	,	"N1"	,
"PbeMT16"	,	"N1"	,
"PbeMT17"	,	"N1"	,
"PbeMT18"	,	"N1"	,
"PbeMT19"	,	"N1"	,
"PbeMT20"	,	"N1"	,
"PbeMT21"	,	"N1"	,
"PbeMT22"	,	"N2"	,
"PbeMT23"	,	"N2"	,
"PbeMT24"	,	"N2"	,
"PbeMT25"	,	"N2"	,
"PbeMT26"	,	"N2"	,
"PbeMT27"	,	"N2"	,
"PbeMT28"	,	"N3"	,
"PbeMT29"	,	"N3"	,
"PbeMT30"	,	"N3"	,
"PbeMT31"	,	"N3"	,
"PbeMT32"	,	"N3"	,
"PbeMT33"	,	"N3"	,
"PbeMT34"	,	"N3"	,
"PbeMT35"	,	"N3"	,
"PbeMT36"	,	"N3"	,
"PbeMT37"	,	"N3"	,
"PbeMT38"	,	"N4"	,
"PbeMT39"	,	"N4"	,
"PbeMT40"	,	"N4"	,
"PbeMT41"	,	"N4"	,
"PbeMT42"	,	"N4"	,
"PbeMT43"	,	"N4"	,
"PbeMT44"	,	"N4"	,
"PbeMT45"	,	"N5"	,
"PbeMT46"	,	"Sunda"	,
"PbeMT47"	,	"Sunda"	,
"PbeMT48"	,	"Sunda"	,
"PbeMT49"	,	"Sunda"	,
"PbeMT50"	,	"Sunda"	,
"PbeMT51"	,	"Sunda"	,
"PbeMT52"	,	"Sunda"	,
"PbeMT53"	,	"Sunda"	,
"PbeMT54"	,	"Sunda"	,
"PbeMT55"	,	"Sunda"	,
"PbeMT56"	,	"N1"	,
"PbeMT57"	,	"N4"	,
"PbeMT58"	,	"N4"	,
"PbeMT59"	,	"N4"	,
"MtHapCode"	,	"MtHapGroup"
);

while (<I>) {
	s/\r\n//;
	my @a = split /\t/;
	$a[5] = $group{$a[4]};
	print O join("\t", @a), "\n";
}
close I;
close O;
