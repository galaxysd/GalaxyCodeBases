#!/usr/bin/perl
use strict;
use warnings;

open I75, "<", "scaffold188SNP.16.gt";
open O, ">", "scaffold188_linkage_map.svg";

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
  foreach (13..19) {
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

print O "<?xml version=\"1.0\"?>\n\n<svg xmlns=\"http://www.w3.org/2000/svg\">\n\n";

my $yy;
foreach my $a (@ra75) {
  ++$yy;
  my $y = 10*$yy;
  print O "<text x=\"0\" y=\"$y\" font-family=\"Courier\" font-size=\"8\" fill=\"black\">${$a}[0]</text>\n";
  print O "<text x=\"70\" y=\"$y\" font-family=\"Courier\" font-size=\"8\" fill=\"black\">${$a}[1]</text>\n";
  foreach (2..17) {
    my $xt = 30*$_ + 60;
    my $yt = 10*$yy;
    my $xr = 30*$_ + 50;
    my $yr = 10*$yy - 8;
    my $h = length ${$a}[$_];
    my $color;
    if ($h == 2) {
      $color = "green";
    } elsif ($h == 1) {
      if (${$a}[18] eq "equal") {
	$color = "pink";
      } elsif (${$a}[18] eq ${$a}[$_]) {
	$color = "yellow";
      } else {
	$color = "blue";
      }
    } else {
      $color = "white";
    }
    print O "<rect x=\"$xr\" y=\"$yr\" width=\"30\" height=\"10\" fill=\"$color\" stroke-width=\"0.5\" stroke=\"black\"/>\n";
    print O "<text x=\"$xt\" y=\"$yt\" font-family=\"Courier\" font-size=\"8\" fill=\"black\">${$a}[$_]</text>\n";    
  }
}

print O "\n</svg>\n";
close O;
