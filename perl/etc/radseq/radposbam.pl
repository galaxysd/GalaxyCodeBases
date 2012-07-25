#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
#use Data::Dump qw(ddx);
use Galaxy::IO;
use Galaxy::IO::FASTAQ;

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
my @aux = undef;
my ($name, $comment, $seq, %RefSeq, $ret);
my ($n, $slen) = (0, 0);
while ($ret = &readfq($FHref, \@aux)) {
	($name, $comment, $seq) = @$ret;
	++$n;
	$slen += length($seq);
	$RefSeq{$name} = $seq;
	#warn "$name, $comment, $seq\n";
}
warn "Ref: $n seq of $slen bp.\n";

my %bamFH;
for my $name (@bamfs) {
	$bamFH{$name} = openpipe('samtools view -f64 -F1796',$name);	# +0x40 -0x704
}

