#!/usr/bin/perl
use strict;
use warnings;

open I, "-|", "bcftools view tigers.bcgv.bcf";
open NOTDAD, ">", "notDAD.vcf";
open NOTMUM, ">", "notMUM.vcf";
open LOG, ">", "parenthood.log";

my ($yDad, $nDad, $yMum, $nMum) = (0, 0, 0, 0);
while (<I>) {
  next if /^#/;
  my @line = split /\t/;
  next if $line[7] =~ /^INDEL/;
  next if $line[5] < 999;
  my @JHH001 = split /:/, $line[9];
  next if $JHH001[4] < 99;
  my @BXH011 = split /:/, $line[18];
  next if $BXH011[4] < 99;
  my @GZXJ36 = split /:/, $line[25];
  next if $GZXJ36[4] < 99;
  my @aJ = split /\//, $JHH001[0];
  my @aB = split /\//, $BXH011[0];
  my @aG = split /\//, $GZXJ36[0];
  my ($GB, $GJ);
  foreach my $nG (0..1) {
    foreach my $nB (0..1) {
      ++$GB if $aG[$nG] eq $aB[$nB];
    }
    foreach my $nJ (0..1) {
      ++$GJ if $aG[$nG] eq $aJ[$nJ];
    }
  }
  if ($GB) {
    ++$yDad;
  } else {
    ++$nDad;
    foreach (0..8) {
      print NOTDAD $line[$_], "\t";
    }
    print NOTDAD $line[9], "\t", $line[18], "\t", $line[25], "\n";
  }
  if ($GJ) {
    ++$yMum;
  } else {
    ++$nMum;
    foreach (0..8) {
      print NOTMUM $line[$_], "\t";
    }
    print NOTMUM $line[9], "\t", $line[18], "\t", $line[25], "\n";
  }
}

close I;
close NOTDAD;
close NOTMUM;

print LOG "yBXH011\tnBXH011\tyJHH001\tnJHH001\n";
print LOG "$yDad\t$nDad\t$yMum\t$nMum\n";

close LOG;
