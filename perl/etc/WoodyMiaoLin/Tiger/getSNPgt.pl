#!/usr/bin/perl
use strict;
use warnings;

my $s = shift;
open I, "-|", "bcftools view -I tigers.bcgv.bcf";
open O1, ">", "${s}SNP.18.vcf";
open O2, ">", "${s}SNP.18.gt";

while (<I>) {
  next unless /^$s\t/;
  chomp;
  my @line = split /\t/;
  next if $line[5] < 20;
  $line[7] =~ /;FQ=([0-9\-.]+)/;
  next unless $1 > 0;
  my $a = 0;
  my (@a, @b);
  foreach (9..26) {
    my @sample = split /:/, $line[$_];
    if ($sample[4] >= 20 and $sample[2] >= 1) {
      ++$a;
      push @a, $sample[0];
      push @b, $sample[0];
    } else {
      push @b, "N/A";
      next;
    }
  }
  next if $a < 18;
  next if (!grep {$_ ne $a[0]} @a);
  print O1 "$_\n";
  print O2 "$line[0]\t$line[1]\t$line[3]\t$line[4]\t", (join "\t", @b), "\n";
}

close I;
close O1;
close O2;
