#!/usr/bin/perl -w

use strict;


unless(@ARGV){
	print "\\n\tto calculate the chromsome length and effective length\n";
	print "\n\tperl $0 <fasta file> < N size, default 1> <output file>\n\n";
	exit 0;
}


my $Nsize = $ARGV[1] || 1;
my $outputFile = "$ARGV[2]";

my $name = '';
my %fasta = ();

open IN,"$ARGV[0]" or die "$!\n";
while(my $line = <IN>){
	chomp $line;
	if($line =~ />(.*$)/){
		$name = (split /\s+/,$1)[0];
	}else{
		$fasta{$name} .= $line;
	}
}
close IN;


open OUT,">$outputFile" or die "$!\n";
warn "Chromosome\tLen1\tEffectiveLen\tGapCount\tGapSize\n";
foreach my $name(sort keys %fasta){
	
	my $chrLen = length $fasta{$name};
	my $gapSize = 0;
	my $gapCount = 0;
	my $length = 0;
	while($fasta{$name} =~ /N{$Nsize,}/ig){
		$length = $+[0] - $-[0];
		$gapSize += $length;
		my $start = $-[0] + 1;
		print OUT "$name\t$start\t$+[0]\t$length\n";
		$gapCount++;
	}
	my $effectiveLen = $chrLen - $gapSize;
	warn "$name\t$chrLen\t$effectiveLen\t$gapCount\t$gapSize\n";
}
close OUT;
