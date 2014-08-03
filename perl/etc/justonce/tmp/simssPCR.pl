#!/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

my $cycle = 5;
my (%PCR,@Solution);
push @Solution,0;

for my $i ( 1 .. $cycle ) {
	my @t =  @Solution;
	push @Solution,(1+$_) for @t;
}

++$PCR{$_} for @Solution;
ddx \%PCR;
