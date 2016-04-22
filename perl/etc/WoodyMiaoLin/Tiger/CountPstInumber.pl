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


my %thread;
while (<>) {
  chomp;
  $thread{$_} = threads->new (\&SiteCover, $_);
}
my $t = localtime;
print STDERR "$t\tDestribute threads completed!\n";


my %all;
foreach (keys %thread) {
  $all{$_} = $thread{$_}->join;
}
$t = localtime;
print STDERR "$t\tRead all file completed!\n";

my %stat_site;
foreach my $i (keys %all) {
  $i =~ /^(\d\d)/;
  my $bit = 18 - $1;
  foreach my $s (keys %{$all{$i}}) {
    foreach my $c (keys %{$all{$i}{$s}}) {
      ++$stat_site{$s}{$c}[0] if $all{$i}{$s}{$c}[0];
      ++$stat_site{$s}{$c}[1] if $all{$i}{$s}{$c}[1];
    }
  }
}

my $PstI;
foreach my $s (keys %stat_site) {
  foreach my $c (keys  %{$stat_site{$s}}) {
    ++$PstI if ($stat_site{$s}{$c}[0] and $stat_site{$s}{$c}[1]);
  }
}

$t = localtime;
print STDERR "$t\tCount PstI complete!\n";
print STDERR "Number of PstI is $PstI\n";




