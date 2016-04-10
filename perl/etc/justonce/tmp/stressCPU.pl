#!/usr/bin/perl
use strict;
use warnings;
use threads;

my $t = shift;

my %thread;
foreach (1 .. $t) {
	$thread{$_} = threads->new(\&test, $_);
}

my %count;
foreach (sort keys %thread) {
	$count{$_} = $thread{$_}->join;
}

sub test {
	my $in = shift;
	my $a;
	while ($in) {
		++$a;
		--$a;
	}
}
