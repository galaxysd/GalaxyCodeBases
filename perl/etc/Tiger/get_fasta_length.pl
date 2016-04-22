#!/usr/bin/env perl
use strict;
use warnings;

my %length;
while (<>) {
  chomp;
  s/^>//;
  $/=">";
  my $seq = <>;
  chomp $seq;
  $seq =~ s/\s//g;
  $/="\n";
  $length{$_} = length $seq;
}

foreach (sort keys %length) {
  print "$_\t$length{$_}\n";
}
