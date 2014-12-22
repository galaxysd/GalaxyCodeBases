#!/bin/env perl
use strict;
use warnings;
use Compress::Zlib;

die "Usage: $0 <in_bam> [out_prefix]\n" if @ARGV < 1;
my ($in,$out)=@ARGV;
unless ($out) {
	$out='grep_'.$in;
	$out =~ s/\.[sb]am(\.gz)?//g;
}

open( IN,"-|","samtools view -h $in") or die "Error opening $in: $!\n";
my $OUT = gzopen($out.'.sam.gz', "wb") or die "Cannot open $out.sam.gz: $gzerrno\n" ;
while (my $line = <IN>) {
	if ($line =~ /^\@/) {
		$OUT->gzwrite($line) or die "error writing 1: $gzerrno\n" ;
	} else {
		#my ($id, $flag, $ref, $pos, $mapq, $CIAGR, $mref, $mpos, $isize, $seq, $qual, @OPT) = split /\t/,$line;
		#print "$id, $flag, $ref, $pos, $mapq, $CIAGR, $mref, $mpos, $isize\n";
		my @Dat1 = split /\t/,$line;
		my $line2 = <IN>;
		die '[x]SAM/BAM file not paired !' unless defined($line2);
		my @Dat2 = split /\t/,$line2;
		#if ( $ref eq 'chrEBV' or ($flag & 12) ) {
		if ( $Dat1[2] eq 'chrEBV' or ($Dat1[1] & 12) or $Dat2[2] eq 'chrEBV' or ($Dat2[1] & 12) ) {
			if ($Dat1[0] ne $Dat2[0]) {
				die "[x]SAM/BAM file not paired as $Dat1[0] != $Dat2[0] !\n";
			}
			$OUT->gzwrite($line) or die "error writing 2: $gzerrno\n" ;
			$OUT->gzwrite($line2) or die "error writing 2: $gzerrno\n" ;
		} elsif ( $Dat1[2] eq '*' or $Dat2[2] eq '*' ) {
			$OUT->gzwrite($line) or die "error writing 3: $gzerrno\n" ;
			$OUT->gzwrite($line2) or die "error writing 3: $gzerrno\n" ;
			warn "-->[$line$line2";
		}
	}
}
$OUT->gzclose;
close IN;
