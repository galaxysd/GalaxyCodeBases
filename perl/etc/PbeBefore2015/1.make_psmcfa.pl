#!/usr/bin/perl
use strict;
use warnings;

open I, "-|", "gzip -dc $ARGV[0]";
open O, ">", "$ARGV[1]";

my $mis;
my $ti;
my $tv;
my $dot;
my $one;
my $o;
while (<I>) {
	last unless /(chr\w\d)/;
	print O ">$1\n";
	$/ = ">";
	my $seq = <I>;
	chomp $seq;
	$seq =~ s/\s//g;
	$/ = "\n";
	my $len = length $seq;
	my $t;
	for (my $i=0; $i < $len; $i += 100) {
		my $u = substr $seq, $i, 100;
		last if length($u) < 100;
		my %v = ('N',0,'R',0,'Y',0,'M',0,'K',0,'S',0,'W',0);
		foreach my $j (0 .. 99) {
			++$v{substr($u, $j, 1)};
		}
		my @v = keys %v;
		die "$ARGV[0]\t@v" if @v > 11;
		$mis += $v{N};
		my $ti100 = $v{R} + $v{Y};
		my $tv100 = $v{M} + $v{K} + $v{S} + $v{W};
		$ti += $ti100;
		$tv += $tv100;
		if ($v{"N"} >= 90) {
			print O ".";
			++$dot;
		} elsif ($ti100 or $tv100) {
			print O "1";
			++$one;
		} else {
			print O "0";
			++$o;
		}
		++$t;
		if ($t == 100) {
			print O "\n";
			$t = 0;
		}
	}
}
close I;
close O;
warn "$ARGV[0]\t$mis\t$ti\t$tv\t$ARGV[1]\t$dot\t$one\t$o\n";
