#!/bin/env perl
use strict;
use warnings;

die "Usage: $0 <sorted_bam> <outprefix>\n" if @ARGV != 2;
my ($in,$out)=@ARGV;

open( HEAD,"-|","samtools view -H $in") or die "Error opening $in: $!\n";
my ($LenT,@Heads,@ChrID,%ChrLen,%ChrFlag,%ChrStat);
while (<HEAD>) {
	chomp;
	if (/\@SQ\tSN:(\S+)\tLN:(\d+)/) {
		push @ChrID,$1;
		$LenT += $2;
		$ChrLen{$1} = $2;
		$ChrFlag{$1} = 0;
		$ChrStat{$1} = 0;
		print ">$1: $2\n";
	} else {
		push @Heads,$_;
	}
}
$ChrLen{'*'} = -1;
close HEAD;

open LOG,'>',"$out.log" or die "Error opening $out.log: $!\n";
sub initOUT($) {
	my $outChr = $_[0];
	open( OUT,"|-","samtools view -bS - > $out.$outChr.bam") or die "Error opening $out.$outChr.bam: $!\n";
	print OUT "\@SQ\tSN:$outChr\tLN:$ChrLen{$outChr}\n",join("\n",@Heads),"\n";
}

my $outChr = $ChrID[0];
$ChrFlag{$outChr} = 1;
initOUT($outChr);

open( IN,"-|","samtools view $in") or die "Error opening $in: $!\n";
while (<IN>) {
	my $chr=(split /\t/)[2];
	++$ChrStat{"\t"};
	if ($chr eq '*') {
		++$ChrStat{$chr};
		next;	# unmap are skipped outputing
	}
	if ($chr eq $outChr) {
		print OUT $_;
	} else {
		$ChrFlag{$outChr} = 1;
		close OUT;
		die "File not sorted:[$chr] was before [$outChr].\n" if $ChrFlag{$chr};
		initOUT($chr);
		$outChr = $chr;
		print OUT $_;
	}
	++$ChrStat{$chr};
}
close IN;
close OUT;
print LOG "ChrID\tLen\tLenRatio\tCount\tRatio\tReletiveRatio\n";
for my $chr (@ChrID,'*') {
	print LOG join("\t",$chr,$ChrLen{$chr},$ChrLen{$chr}/$LenT,$ChrStat{$chr},$ChrStat{$chr}/$ChrStat{"\t"},($ChrStat{$chr}/$ChrStat{"\t"})/($ChrLen{$chr}/$LenT)),"\n";
}
close LOG;
