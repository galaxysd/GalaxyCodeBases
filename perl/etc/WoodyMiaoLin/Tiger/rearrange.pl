#!/usr/bin/perl
use strict;
use warnings;

my @ra;
while (<>) {
  chomp;
  my @a = split /\t/;
  my @b;
  $b[0] = $a[0];
  $b[1] = $a[1];
  foreach (4..21) {
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
  push @ra, \@b;
}

foreach (@ra) {
  my $a = join "\t", @{$_};
  print $a, "\n";
}
  
print O "\n</svg>\n";
close O;
