#!/usr/bin/perl
use strict;
use warnings;

my %align;
while (<>) {
  my @a = split /\t/;
  $align{$a[1]}{$a[6]}{$a[8]} = $_;
}

foreach my $c (sort keys %align) {
  foreach my $d (sort {$a <=> $b} keys %{$align{$c}}) {
    foreach my $e (sort {$a <=> $b} keys %{$align{$c}{$d}}) {
      print $align{$c}{$d}{$e};
    }
  }
}
