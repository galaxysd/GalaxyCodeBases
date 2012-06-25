#!/bin/env perl
use strict;
use warnings;
use Galaxy::IO::FASTAQ;
use Galaxy::IO;
use Galaxy::Casava::Eamss qw(doEamss);

die "Usage: $0 <ori-fq> <B-masked-fq> [outprefix].(log|stat)\n" if @ARGV <2;
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
my ($maskedQ,$mlen,$rlen,$t,$ret,%Count,%MaskLen,%Readlen,%ReadType);
my ($totBase,$totMaskedBP,$totMaskedReads)=(0,0,0);
my ($totBaseC,$totMaskedBPC,$totMaskedReadsC)=(0,0,0);
while (1) {
	$dat1 = readfq($fha, \@aux1);	# [$name, $comment, $seq, $qual]
	$dat2 = readfq($fhb, \@aux2);
	if ($dat1 && $dat2) {
		($maskedQ,$mlen) = @{doEamss($$dat1[2],$$dat1[3])};
		if ($$dat1[1] =~ /(\d):([YN]):\d+:([ATCG]+)?/) {	# 1:N:0:TGACCA
			$t = ($1 eq '2')? 1 : 0;
			$t |= 2 if ($2 eq 'Y');
			$t |= 4 if $3;
			++$ReadType{$t};	# 1:Read2, 2:Casava-filtered, 3:with-index
		}
		#print "3-",join('|',$$dat1[0],$$dat1[1],$mlen,$t),"\n$$dat1[2]\n$$dat1[3]\n$maskedQ\n";
		$ret = doCMP($dat1,$dat2,$maskedQ);
		++$Count{$ret};
		++$MaskLen{$mlen};
		$rlen = length $$dat1[2];
		++$Readlen{$rlen};
		++$CountPairs;
		$totBase += $rlen;
		if ($mlen) {
			++$totMaskedReads;
			$totMaskedBP += $mlen;
			unless ($t | 2) {
				$totMaskedBPC += $mlen;
				++$totMaskedReadsC;
			}
		}
		$totBaseC += $rlen unless ($t & 2);
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
open STAT,'>',"${out}.stat" or die "Error opening $out.stat:$!\n";
my $str = "Out Pairs: $CountPairs\nFQ1 over hang: $Count1\nFQ2 over hang: $Count2\n";
print $str;
print LOG "From [$inA]&[$inB] to [$out.stat]\n$str";

$totMaskedReads = -1 unless $totMaskedReads;
$t = int(0.5+10*$totMaskedBP/$totMaskedReads)/10;
$str = "\nTotal Reads, Base: $CountPairs, $totBase\tMasked Reads, Base: $totMaskedBP, $totMaskedReads\tAvg Masked Length: $t\n";
$totMaskedReadsC = -1 unless $totMaskedReadsC;
$t = int(0.5+10*$totMaskedBPC/$totMaskedReadsC)/10;
$str .= "After Casava-filter:\nTotal Base: $totBaseC\tMasked Reads, Base: $totMaskedBPC, $totMaskedReadsC\tAvg Masked Length: $t\n\n";
print $str;
print LOG $str,"@ $totBase\t$totMaskedBP\t$inA\n";

$str="[Compared] (0 or 8 is OK)\n";
for my $k (sort {$a <=> $b} keys %Count) {
	$str .= sprintf("%#06b\t(%#x)\t%d\n",$k,$k,$Count{$k});
}
print $str;
print STAT $str;

$str="\n[Read_Type] (1:Read2, 2:Casava-filtered, 4:with-index)\n";
for my $k (sort {$a <=> $b} keys %ReadType) {
	$str .= sprintf("%#06b\t(%#x)\t%d\n",$k,$k,$ReadType{$k});
}
print $str;
print STAT $str;

$str="\n[Read_Length]\n";
for my $k (sort {$a <=> $b} keys %Readlen) {
	$str .= sprintf("%d\t%d\n",$k,$Readlen{$k});
}
#print $str;
print STAT $str;

$str="\n[Masked_Length]\n";
for my $k (sort {$a <=> $b} keys %MaskLen) {
	$str .= sprintf("%d\t%d\n",$k,$MaskLen{$k});
}
#print $str;
print STAT $str;

close LOG;
close STAT;


sub doCMP($$$) {
	my ($dat1,$dat2,$maskedQ) = @_;
	my $flag = 0;
	for my $i (0 .. 3) {
		$flag |= 1<<$i if $$dat1[$i] ne $$dat2[$i];
	}
	$flag |= 1<<4 if $maskedQ ne $$dat2[3];
	return $flag;
}