#!/usr/bin/perl
use strict;
use warnings;

open I, "<", "scaffold188SNP.18.gt";
open O, ">", "scaffold188SNP.16.gt";

while (<I>) {
  my @a = split /\t/;
  splice @a, 19, 2;
  print O join "\t", @a;
}

close I;
close O;

