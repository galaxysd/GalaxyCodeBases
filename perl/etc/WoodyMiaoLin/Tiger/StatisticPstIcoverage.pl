#!/usr/bin/perl
use strict;
use warnings;
use threads;


sub SiteCover {
  my $file = shift;
  my %count;
  open my $in, "-|", "samtools view -f 64 -F 1796 $file" or die "Error opening $file: $!\n";
  while (<$in>) {
    my @align = split /\t/;
    my @match = $align[5] =~ /(\d+)M/g;
    my $sum_match = 0;
    $sum_match += $_ foreach @match;
    next if $sum_match < 64;
    @match = sort {$b <=> $a} @match;
    next if $match[0] < 30;
    if ($align[1] & 0x10) {
      next if $align[5] =~ /[NHP=X]/;
      next unless $align[5] =~ /(\d+)M$/;
      next unless $1 >= 5;
      my @del = $align[5] =~ /(\d+)D/g;
      my $sum_del = 0;
      $sum_del += $_ foreach @del;
      my $c = $align[3] + $sum_match + $sum_del - 1;
      $count{$align[2]}{$c}[0] = 1;
    } else {
      next unless $align[5] =~ /^(\d+)M/;
      next unless $1 >= 5;
      $count{$align[2]}{$align[3]+3}[1] = 1;
    }
  }
  close $in;
  my $t = localtime;
  print STDERR "$t\tCount file $file complete!\n";
  return \%count;
}


my %thread; #Creat a hash of threads.
while (<>) {
  chomp;
  $thread{$_} = threads->new (\&SiteCover, $_);
}
my $t = localtime;
print STDERR "$t\tDestribute threads completed!\n";


my %all; #Threads return a hash of all individual's coverage.
foreach (keys %thread) {
  $all{$_} = $thread{$_}->join;
}
$t = localtime;
print STDERR "$t\tRead all file completed!\n";


my (%stat_site, %stat_indiv); #Count sites' coverage and individuals' coverage.
foreach my $i (keys %all) {
  $i =~ /^(\d\d)/;
  my $bit = 18 - $1;
  foreach my $s (keys %{$all{$i}}) {
    foreach my $c (keys %{$all{$i}{$s}}) {
      $stat_site{$s}{$c}[0] = 0 unless defined $stat_site{$s}{$c}[0];
      $stat_site{$s}{$c}[1] = 0 unless defined $stat_site{$s}{$c}[1];
      if (defined $all{$i}{$s}{$c}[0]) {
	++$stat_indiv{$i}[0];
	$stat_site{$s}{$c}[0] = $stat_site{$s}{$c}[0] | ($all{$i}{$s}{$c}[0] <<= $bit);
      } else {
	$all{$i}{$s}{$c}[0] = 0;
      }
      if (defined $all{$i}{$s}{$c}[1]) {
	++$stat_indiv{$i}[1];
        $stat_site{$s}{$c}[1] = $stat_site{$s}{$c}[1] | ($all{$i}{$s}{$c}[1] <<= $bit);
      } else {
        $all{$i}{$s}{$c}[1] = 0;
      }
      $all{$i}{$s}{$c}[2] = $all{$i}{$s}{$c}[0] | $all{$i}{$s}{$c}[1];
      $all{$i}{$s}{$c}[3] = $all{$i}{$s}{$c}[0] & $all{$i}{$s}{$c}[1];
      ++$stat_indiv{$i}[2] if $all{$i}{$s}{$c}[2];
      ++$stat_indiv{$i}[3] if $all{$i}{$s}{$c}[3];
    }
  }
}
$t = localtime;
print STDERR "$t\tDo stastistics completed!\n";


my @stat_time; #Print site coverage and count coverage of times.
open Osite, ">", "PstIcoverage.lst";
print Osite "#SCAFFOLD\tCOORDINATE\tREVERSE\tFORWARD\tUNION\tINTERSECTION\n";
foreach my $s (sort keys %stat_site) {
  foreach my $c (sort {$a <=> $b} keys %{$stat_site{$s}}) {
    $stat_site{$s}{$c}[2] = $stat_site{$s}{$c}[1] | $stat_site{$s}{$c}[0];
    $stat_site{$s}{$c}[3] = $stat_site{$s}{$c}[1] & $stat_site{$s}{$c}[0];
    printf Osite "%s\t%d\t%018b\t%018b\t%018b\t%018b\n", $s, $c, @{$stat_site{$s}{$c}};
    foreach my $r (0..3) {
      my @n = split //, sprintf("%b", $stat_site{$s}{$c}[$r]);
      my $n = 0;
      $n += $_ foreach @n;
      ++$stat_time[$n][$r];
    }
  }
}
close Osite;
$t = localtime;
print STDERR "$t\tPrint to PstIcoverage.lst completed!\n";


open Oindiv, ">", "PstIindivCoverage.lst"; #Print individual and time coverage.
print Oindiv "#INDIV.\tREVERSE\tFORWARD\tUNION\tINTERSECTION\n";
foreach my $i (sort keys %stat_indiv) {
  printf Oindiv "%s\t%d\t%d\t%d\t%d\n", $i, @{$stat_indiv{$i}};
}
print Oindiv "#TIME\tREVERSE\tFORWARD\tUNION\tINTERSECTION\tTIME\tREVERSE\tFORWARD\tUNION\tINTERSECTION\n";
my @accum;
foreach my $n1 (1..18) {
  foreach my $n2 (0..3) {
    foreach my $n3 ($n1..18) {
      $accum[$n1][$n2] += $stat_time[$n3][$n2];
    }
  }
  printf Oindiv "%s\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\n", ">=$n1", @{$accum[$n1]}, "=$n1", @{$stat_time[$n1]};
}
close Oindiv;
$t = localtime;
print STDERR "$t\tPrint to PstIindivCoverage.lst completed!\n";


$t = localtime;
print STDERR "$t\tDone!\n";


