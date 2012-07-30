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

die "Usage: $0 <reference> <out> <sorted bam files>\n" if @ARGV<2;
my $fa=shift;
my $out=shift;
my @bamfs = @ARGV;

my $ReadLen = 101-5;
my $minOKalnLen = 30;
my $minAlignLen = int($ReadLen * 0.6);
$minAlignLen = $minOKalnLen if $minAlignLen < $minOKalnLen;
die if $ReadLen < $minOKalnLen;

my $Eseq="CTGCAG";
my $EcutAt=5;
my $EseqL="CTGCA";#substr $Eseq,0,$EcutAt;
my $EseqR="TGCAG";

my $EfwdTerm5=1-$EcutAt+1;	# -3
my $EfwdTerm3=1-$EcutAt+$ReadLen;	# 91
my $ErevTerm3=0;	# 0
my $ErevTerm5=$ErevTerm3-$ReadLen+1;	# -94
# 12008 -> [12005,12099],[11914,12008]
my $Rfwd2Ec = -$EfwdTerm5;	# 3
my $Rrev2Ec = -$ErevTerm3;	# 0

my $t;
open O,'>',$out or die "Error opening $out : $!\n";
$t = "# Ref: [$fa], Enzyme: [$Eseq], Cut after $EcutAt\n# Bams: [@bamfs]\n\n";
print O $t;
print $t;

my $FHref = openfile($fa);
my @aux = undef;
my ($name, $comment, $seq, %RefSeq, %RefLen, $ret);
my ($n, $slen) = (0, 0);
while ($ret = &readfq($FHref, \@aux)) {
	($name, $comment, $seq) = @$ret;
	++$n;
	$RefSeq{$name} = $seq;
	$RefLen{$name} = length($seq);
	$slen += $RefLen{$name};
	#warn "$name, $comment, $seq\n";
}
warn "Ref: $n seq of $slen bp.\n";
my $th = openpipe('samtools view -H',$bamfs[0]);
my @RefOrder;
while (<$th>) {
	/^\@SQ\tSN:([^ ]+) LN:(\d+)$/ or next;
	push @RefOrder,$1;
	die "[x]Chr:$1 Len:$2 != $RefLen{$name}.\n" if $2 != $RefLen{$name};
}

sub cigar2reflen($) {
	my $cigar = $_[0];
	my @cigar = $cigar =~ /(\d+)(\w)/g;
	my ($reflen,$maxM)=(0,0);
	while (@cigar) {
		my ($len,$op) = splice(@cigar,0,2);
		if ($op eq 'M') {
			$reflen += $len;
			$maxM = $len if $maxM < $len;
			next;
		}
		$reflen += $len if $op eq 'D';
	}
	return [$reflen,$maxM];
}

my (@bamData,@ta);
my $fhflag = 0;
for my $name (@bamfs) {
	$th = openpipe('samtools view -f64 -F1796',$name);	# +0x40 -0x704
	my $t = readline $th;
	#my ($id,$flag,$ref,$pos,$mapq,$cigar,$MRNM,$MPOS,$ISIZE,$seq,$qual,@opt) = split /\t/,$t;
	@ta = split /\t/,$t;
	push @bamData,[$th,1,$ta[2],$ta[3],\@ta];
	++$fhflag;
}

while ($fhflag > 0) {
	@bamData = sort { $a->[2] <=> $b->[2] } @bamData;
}
