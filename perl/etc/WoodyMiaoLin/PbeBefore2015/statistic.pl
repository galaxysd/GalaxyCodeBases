#!/usr/bin/perl
use strict;
use warnings;

$/ = "\r\n";
$_ = <>;
chomp;
my @locus_id = split /\t\t/;
$locus_id[0] =~ s/Animal\tLocale of origin\tRegion of origin\t//;

my (%sample, @locus);
while (<>) {
  chomp;
  my @a = split /\t/;
  my $n;
  $sample{$a[0]} = 0;
  for (my $i = 4; $i <= 36; $i += 2) {
    if ($a[$i-1]) {
      ++$sample{$a[0]};
      ++$locus[0.5*$i];
    }
  }
}

print "#Sample\tNo.loci\n";
foreach (sort {$sample{$b} <=> $sample{$a}} keys %sample) {
  print "$_\t$sample{$_}\n";
}

print "#Locus\tNo.sample\n";
foreach (0..16) {
  print "$locus_id[$_]\t$locus[$_+2]\n"
}
