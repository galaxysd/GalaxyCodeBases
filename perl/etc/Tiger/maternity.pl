#!/usr/bin/perl
use strict;
use warnings;

open I, "-|", "bcftools view tigers.bcgv.bcf";

my (@isMom, @noMom);
while (<I>) {
  next if /^#/;
  my @line = split /\t/;
  next if $line[7] =~ /^INDEL/;
  next if $line[5] < 20;
  my @JHHgt = split /:/, $line[9];
  next if $JHHgt[4] < 20;
  my ($Ja1, $Ja2) = split /\//, $JHHgt[0];
  foreach (10..26) {
    my @ANYgt = split /:/, $line[$_];
    next if $ANYgt[4] < 20;
    my ($Aa1, $Aa2) = split /\//, $ANYgt[0];
    my $a;
    ++$a if $Aa1 eq $Ja1;
    ++$a if $Aa1 eq $Ja2;
    ++$a if $Aa2 eq $Ja1;
    ++$a if $Aa2 eq $Ja2;
    $_ -= 10;
    if ($a) {
      ++$isMom[$_];
    } else {
      ++$noMom[$_];
    }
  }
}

close I;
open O, ">", "maternity.log";

print O "#\tGZXJ03\tGZXJ05\tGZXJ26\tGZXJ27\tGZXJ28\tGZXJ29\tGZXJ30\tGZXJ33\tBHX011\tBHX019\tGZXJ04\tGZXJ06\tGZXJ31\tGZXJ32\tGZXJ36\tGZXJ37\tGZXJ38\n";
print O "isMom\t";
print O join "\t", @isMom;
print O "\nnoMom\t";
print O join "\t", @noMom;
print O "\n";

close O;
