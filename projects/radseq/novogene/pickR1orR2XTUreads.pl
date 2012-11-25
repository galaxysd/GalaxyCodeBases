#!/usr/bin/env perl
use strict;
use warnings;

die "Usage: $0 <in sam.gz> <out sam.gz>\n" if @ARGV<2;

my $in=shift;
my $out=shift;

open I, '-|', "/bin/gzip -dc $in" or die "Error opening $in: $!\n";
open O, '|-', "/bin/gzip -9c >$out" or die "Error opening $out: $!\n";
open LOG, '>>', 'pickR1orR2XTUreads.log';

while (<I>) {
  if (/^@/) {
    print O;
    next;
  }
  my $r1 = $_;
  my @r1 = split /\t/, $r1;
  my $r2 = <I>;
  my @r2 = split /\t/, $r2;
  if (($r1[1] & 0x40) and ($r2[1] & 0x80) and ($r1[0] eq $r2[0])) {
    print O $r1,$r2 if ((grep /XT:A:U/, @r1) or (grep /XT:A:U/, @r2));
  } else {
    print LOG "Error!\t$in\t$r1[0]\t$r1[1]\t$r2[0]\t$r2[1]\n";
    exit;
  }
}

print LOG "Complete!\t$in\t$out\n";

close I;
close O;
close LOG;

