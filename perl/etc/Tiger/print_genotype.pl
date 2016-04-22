#!/usr/bin/perl
use strict;
use warnings;

open I75, "<", "scaffold188SNP.18.gt";
open O, ">", "scaffold188.18.genotype";

my @ra75;
while (<I75>) {
  chomp;
  my @a = split /\t/;
  my @b;
  $b[0] = $a[0];
  $b[1] = $a[1];
  foreach (4..12) {
    my ($a1, $a2) = split /\//, $a[$_];
    if (($a1 eq "N") and ($a2 eq "A")) {
      push @b, "none";
    } elsif (($a1 eq "0") and ($a2 eq "0")) {
      push @b, $a[2];
    } elsif (($a1 eq "1") and ($a2 eq "1")) {
      push @b, $a[3];
    } else {
      push @b, "$a[2]$a[3]";
    }
  }
  my $b0 = 0;
  my $b1 = 0;
  foreach (13..21) {
    my ($a1, $a2) = split /\//, $a[$_];
    if (($a1 eq "N") and ($a2 eq "A")) {
      push @b, "none";
    } elsif (($a1 eq "0") and ($a2 eq "0")) {
      ++$b0;
      push @b, $a[2];
    } elsif (($a1 eq "1") and ($a2 eq "1")) {
      ++$b1;
      push @b, $a[3];
    } else {
      push @b, "$a[2]$a[3]";
    }
  }
  if ($b0 > $b1) {
    push @b, $a[2];
  } elsif ($b0 < $b1) {
    push @b, $a[3];
  } else {
    push @b, "equal"
  }
  push @ra75, \@b;
}
close I75;

foreach (@ra75) {
  my $a = join "\t", @{$_};
  print O $a, "\n";
}

close O;
