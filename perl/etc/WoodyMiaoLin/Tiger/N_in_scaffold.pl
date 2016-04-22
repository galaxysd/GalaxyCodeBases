#!/usr/bin/perl
use strict;
use warnings;

<>;
$/ = ">";
my $_ = <>;
s/\s//g;
my @a = split //;
my $n;
foreach (@a) {
  ++$n;
  if (!/N/) {
    print "$n\n";
  }
}
