#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120627
=cut
use strict;
use warnings;
#use Data::Dump qw(ddx);
use Galaxy::IO::FASTAQ;
use Galaxy::IO;

die "Usage: $0 <tag list> <fq1> <fq2> <out prefix>\n" if @ARGV<4;
my $tagf=shift;
my $fq1f=shift;
my $fq2f=shift;
my $outp=shift;
my $maxmismark = 1;

my (%BarSeq2idn,@BarSeq,@BarSeqs);
open L,'<',$tagf or die "Error opening $tagf:$!\n";
while (<L>) {
	chomp;
	my ($seq,$id,$name)=split /\t/;
	next unless $name;
	$seq = uc $seq;
	push @BarSeq,$seq;
	push @BarSeqs,[split //,$seq];
	$BarSeq2idn{$seq} = [$id,$name,0,0];
}
$BarSeq2idn{'N'} = ['NA','Unknown',0,0];
my $BARLEN = length $BarSeq[0];
#ddx \@BarSeq;
#ddx \%BarSeq2idn;
die if $BARLEN != length $BarSeq[-1];
for my $k (keys %BarSeq2idn) {
	my $fname = join('.',$outp,$BarSeq2idn{$k}->[0],$BarSeq2idn{$k}->[1],'1.fq.gz');
	my ($fha,$fhb);
	open $fha,'|-',"gzip -9c >$fname" or die "Error opening output $fname:$!\n";
	$fname = join('.',$outp,$BarSeq2idn{$k}->[0],$BarSeq2idn{$k}->[1],'2.fq.gz');
	open $fhb,'|-',"gzip -9c >$fname" or die "Error opening output $fname:$!\n";
	($BarSeq2idn{$k}->[2],$BarSeq2idn{$k}->[3]) = ($fha,$fhb);
}
#open $UNKNOWN,'|-',"gzip -9c >$outp.NA.Unknown.fq.gz" or die "Error opening output file:$!\n";

sub CmpBarSeq($$) {
	my ($barseq,$maxmark)=@_;
	return [$BarSeq2idn{$barseq},0] if (exists $BarSeq2idn{$barseq});
	my $ret = [$BarSeq2idn{'N'},-1];
	my @seqss = split //,$barseq;
	my ($minmark,$mismark,$i,%mark2i,%marks)=(999999);
	for ($i=0; $i <= $#BarSeq; ++$i) {
		$mismark = 0;
		for (my $j=0; $j <= $#seqss; ++$j) {
			if ($seqss[$j] eq $BarSeqs[$i][$j]) {
				next;
			} elsif ($seqss[$j] eq 'N') {
				$mismark += 1;
			} else {
				$mismark += 100;
			}
		}
		$mark2i{$mismark} = $i;
		++$marks{$mismark};
		$minmark = $mismark if $minmark > $mismark;
	}
	if ($marks{$minmark} == 1 and $minmark <= $maxmark) {
		$ret = [$BarSeq2idn{$BarSeq[$mark2i{$minmark}]},$minmark];
	} else { $ret->[1] = "$minmark.$marks{$minmark}"; }
	return $ret;
}

my $fha=openfile($fq1f);
my $fhb=openfile($fq2f);
my (@aux1,@aux2);
my ($Count1,$Count2,$CountPairs)=(0,0,0);
my ($dat1,$dat2);
while (1) {
	$dat1 = readfq($fha, \@aux1);	# [$name, $comment, $seq, $qual]
	$dat2 = readfq($fhb, \@aux2);
	if ($dat1 && $dat2) {
		my $bar = substr $$dat1[2],0,$BARLEN;
		my ($kret,$themark) = @{CmpBarSeq($bar,$maxmismark)};
		my $fha = $kret->[2];
		my $fhb = $kret->[3];
		print $fha join("\n",
			join(' ',$$dat1[0],$$dat1[1],$kret->[0],$themark,$kret->[1]),
			$$dat1[2],'+',$$dat1[3]),"\n";
		print $fhb join("\n",
			join(' ',$$dat2[0],$$dat2[1],$kret->[0],$themark,$kret->[1]),
			$$dat2[2],'+',$$dat2[3]),"\n";
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

for my $k (keys %BarSeq2idn) {
	close $BarSeq2idn{$k}->[2];
	close $BarSeq2idn{$k}->[3];
}
