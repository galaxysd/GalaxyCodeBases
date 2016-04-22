#!/usr/bin/perl
use strict;
use warnings;
use threads;


my $t = localtime;
print STDERR "$t\tStart!\n";


my (@indiv, @stat_time, %nSNP, %iSNP);
while (<>) {
  next if /^#/;
  chomp;
  my @line = split /\t/;
  next if $line[7] =~ /^INDEL/;
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
      push @b, $_;
    } else {
      next;
    }
  }
  next if $a < 6;
  next if (!grep {$_ ne $a[0]} @a);
  print "$_\n";
  ++$stat_time[$a];
  $nSNP{$line[0]}{$line[1]} = $a;
  $iSNP{$line[0]}{$line[1]} = 0;
  foreach (@b) {
    ++$indiv[$_];
    $iSNP{$line[0]}{$line[1]} |= (1 << (26-$_));
  }
}
$t = localtime;
print STDERR "$t\tRead file completed!\n";


open OSNP, ">", "SNPcoverage.lst";
print OSNP "#SCAFFOLD\tCOORDINATE\tINVID.\tNO.\n";
foreach my $c (sort keys %nSNP) {
  foreach my $d (sort {$a <=> $b} keys %{$nSNP{$c}}) {
    printf OSNP "%s\t%d\t%018b\t%d\n", $c, $d, $iSNP{$c}{$d}, $nSNP{$c}{$d};
  }
}
close OSNP;
$t = localtime;
print STDERR "$t\tPrint to SNPCoverage.lst completed!\n";


my @id = qw/01_JHH001 02_GZXJ03 03_GZXJ05 04_GZXJ26 05_GZXJ27 06_GZXJ28 07_GZXJ29 08_GZXJ30 09_GZXJ33 10_BHX011 11_BHX019 12_GZXJ04 13_GZXJ06 14_GZXJ31 15_GZXJ32 16_GZXJ36 17_GZXJ37 18_GZXJ38/;
open Oindiv, ">", "SNPindivCoverage.lst";
print Oindiv "#INDIV.\tNO.\n";
foreach (9..26) {
  print Oindiv "$id[$_-9]\t$indiv[$_]\n";
}
print Oindiv "#TIME\tNO.\tTIME\tNO.\n";
my @accum;
foreach my $e (6..18) {
  foreach ($e..18) {
    $accum[$e] += $stat_time[$_];
  }
  printf Oindiv "%s\t%d\t%s\t%d\n", ">=$e", $accum[$e], "=$e", $stat_time[$e];
}
close Oindiv;
$t = localtime;
print STDERR "$t\tPrint to SNPindivCoverage.lst completed!\n";


$t = localtime;
print STDERR "$t\tDone!\n";


  
    
  
