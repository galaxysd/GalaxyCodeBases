#!/bin/env perl
use strict;
use warnings;
use Galaxy::IO::FASTAQ;
use Galaxy::IO;
use Galaxy::Casava::Eamss qw(doEamss);

die "Usage: $0 <ori-fq> <B-masked-fq> <outprefix>.(log|stat)\n" if @ARGV <2;
my ($inA,$inB,$out)=@ARGV;
$out=$inA unless $out;
warn "From [$inA]&[$inB] to [$out].(log|stat)\n";

my $fha=openfile($inA);
my $fhb=openfile($inB);

#my @aux1 = undef;
#my @aux2 = undef;
my (@aux1,@aux2);
my ($Count1,$Count2,$CountPairs)=(0,0,0);

sub doCMP($$$);

my ($dat1,$dat2);
my ($maskedQ,$ret,%Count);
while (1) {
	$dat1 = readfq($fha, \@aux1);	# [$name, $comment, $seq, $qual]
	$dat2 = readfq($fhb, \@aux2);
	if ($dat1 && $dat2) {
		$maskedQ = doEamss($$dat1[2],$$dat1[3]);
		#print "3-",join('|',$$dat1[0],$$dat1[1]),"\n$$dat1[2]\n$$dat1[3]\n$maskedQ\n";
		$ret = doCMP($dat1,$dat2,$maskedQ);
		++$Count{$ret};
		++$CountPairs;
	} elsif ($dat1) {
		++$Count1;
		print "1-",join('|',@$dat1),"\n";
	} elsif ($dat2) {
		++$Count2;
		print "2-",join('|',@$dat2),"\n";
	} else {
		last;
	}
}
close $fha;
close $fhb;

open LOG,'>',"${out}.log" or die "Error opening $out.log:$!\n";
my $str = "Out Pairs: $CountPairs\nFQ1 over hang: $Count1\nFQ2 over hang:$Count2\n";
print $str;
print LOG "From [$inA]&[$inB] to [$out.stat]\n$str";
$str="\nStat (0 or 8 are OK):\n";
print $str;
print LOG $str;
$str='';
for my $k (sort {$a <=> $b} keys %Count) {
	$str .= sprintf("%#06b,%#x\t%d\n",$k,$k,$Count{$k});
}
print $str;
print LOG $str;
close LOG;



sub doCMP($$$) {
	my ($dat1,$dat2,$maskedQ) = @_;
	my $flag = 0;
	for my $i (0 .. 3) {
		$flag |= 1<<$i if $$dat1[$i] ne $$dat2[$i];
	}
	$flag |= 1<<4 if $maskedQ ne $$dat2[3];
	return $flag;
}