#!/usr/bin/perl
use strict;
use warnings;

my $in1 = shift;
my $in2 = shift;
open I1, "<", "$in1";
open I2, "-|", "gzip -dc $in2";

my @notN;
while (<I1>) {
  chomp;
  push @notN, $_;
}
close I1;
#print STDERR "read $in1 complete!\n";

my $n75 = 0;
my $n1458 = 0;
my (@cover, @line, @depth);
while (<I2>) {
  chomp;
  my @a = split /\t/;
  $line[$a[1]] = $_;
  my $c = 0;
  my @d;
  foreach (9..11) {
    my @b = split /:/, $a[$_];
    push @d, $b[1];
    if ($b[1] >= 3) {
      ++$c;
    }
  }
  $depth[$a[1]] = join ",", @d;
  $cover[$a[1]] = 1 if $c == 3;
}
close I2;
#print STDERR "read $in2 complete!\n";

foreach (@notN) {
  $line[$_] = "N/A" unless defined $line[$_];
  $depth[$_] = "0,0,0" unless defined $depth[$_];
  print "$_\t$depth[$_]\t$line[$_]\n" unless defined $cover[$_];
}
