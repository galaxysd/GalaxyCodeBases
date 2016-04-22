#!/usr/bin/perl
use strict;
use warnings;

open I, "<", "Pbe_1_row.txt";
open O, ">", "Pbe_1_column.txt";

my %geo_ind;
while (<I>) {
	chomp;
	my @b = split /\t/;
	my $ind = shift @b;
	my $geo = shift @b;
	$geo_ind{$geo}{$ind} = \@b;
}
close I;

foreach my $b (sort keys %geo_ind) {
	print O "$b\r\n";
	foreach my $c (sort keys %{$geo_ind{$b}}) {
		my @d = @{$geo_ind{$b}{$c}};
		print O "\t\t$c\t1\t$d[0]\t$d[2]\t$d[4]\t$d[6]\t$d[8]\t$d[10]\t$d[12]\t$d[14]\t$d[16]\t$d[18]\t$d[20]\t$d[22]\t$d[24]\t$d[26]\r\n";
		print O "\t\t\t\t$d[1]\t$d[3]\t$d[5]\t$d[7]\t$d[9]\t$d[11]\t$d[13]\t$d[15]\t$d[17]\t$d[19]\t$d[21]\t$d[23]\t$d[25]\t$d[27]\r\n";
	}
}
close O;
