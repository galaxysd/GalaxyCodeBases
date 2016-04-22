#!/usr/bin/perl
use strict;
use warnings;

my $SNP = 0;
while (<>) {
  next if /^#/;
  chomp;
  my @line = split /\t/;
  next if $line[7] =~ /^INDEL/;
  next if $line[5] < 20;
  my $a = 0;
  my @gt;
  foreach (9..11) {
    my @sample = split /:/, $line[$_];
    if ($sample[4] >= 20 and $sample[2] >= 1) {
      ++$a;
      push @gt, $sample[0];
    } else {
      next;
    }
  }
  next if $a < 2;
  my %allele;
  foreach (@gt) {
    my @a =  split /\//;
    $allele{$a[0]} = 1;
    $allele{$a[1]} = 1;
  }
  ++$SNP if (keys %allele) > 1;
}
print "$SNP\n";
