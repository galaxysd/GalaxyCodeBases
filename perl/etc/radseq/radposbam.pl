#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
#use Data::Dump qw(ddx);
use Galaxy::IO;

die "Usage: $0 <reference> <out> <bam files>\n" if @ARGV<2;
my $fa=shift;
my $out=shift;
my @bamfs = @ARGV;

my $Eseq="CTGCAG";
my $EcutAt=5;

my $t;
open O,'>',$out or die "Error opening $out : $!\n";
$t = "# Ref: [$fa], Enzyme: [$Eseq], Cut after $EcutAt\n# Bams: [@bamfs]\n\n";
print O $t;
print $t;

my $FHref = openfile($fa);


