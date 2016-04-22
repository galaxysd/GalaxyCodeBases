#!/usr/bin/perl
use strict;
use warnings;

while (<>) {
  my $a = $_;
  chomp;
  my @a = split /\t/;
  foreach (9..11) {
    my @b = split /:/, $a[$_];
    if ($b[2] != 0) {
      print $a;
      last;
    }
  }
}
