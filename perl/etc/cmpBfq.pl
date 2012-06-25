#!/bin/env perl
use strict;
use warnings;
use Galaxy::IO::FASTAQ;
use Galaxy::IO;

die "Usage: $0 <ori-fq> <B-masked-fq> <outprefix>.(log|stat)\n" if @ARGV <2;
my ($inA,$inB,$out)=@ARGV;
$out=$inA unless $out;
warn "From [$inA]&[$inB] to [$out].(log|stat)\n";

my $fha=openfile($inA);
my $fhb=openfile($inB);

#my @aux1 = undef;
#my @aux2 = undef;
my (@aux1,@aux2);
my ($name, $comment, $seq, $qual);
my ($Count1,$Count2,$CountPairs)=(0,0,0);

my ($dat1,$dat2);
while (1) {
	$dat1 = readfq($fha, \@aux1);
	$dat2 = readfq($fhb, \@aux2);
	if ($dat1 && $dat2) {
		;
		#print "1-",join('|',@$dat1),"\n";
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

#while (($name, $comment, $seq, $qual) = readfq($fha, \@aux)) {

close $fha;
close $fhb;

open LOG,'>',"${out}.log" or die "Error opening $out.log:$!\n";
my $str = "Out Pairs: $CountPairs\nFQ1 over hang: $Count1\nFQ2 over hang:$Count2\n";
print $str;
print LOG "From [$inA]&[$inB] to [$out.fq.xz]\n$str";
close LOG;
