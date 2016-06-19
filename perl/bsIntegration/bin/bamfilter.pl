#!/usr/bin/env perl
use strict;
use warnings;
#use IO::Unread;

die "Usage: $0 <min MapQ> <min softClip> <in.bam> <out.bam>\n" if @ARGV <3;
my ($minMAPQ,$minS,$in,$out) = @ARGV;

open(HEAD,"-|","samtools view -H $in") or die "Error opening $in: $!\n";
open(OUT,"|-","samtools view -b - >$out") or die "Error opening $out: $!\n";
while (<HEAD>) {
	print OUT $_;
}
close HEAD;

open(IN,"-|","samtools view $in") or die "Error opening $in: $!\n";
my ($lastid,@Reads);
while(<IN>) {
	my @Dat = split /\t/;
	if ($lastid and ($lastid ne $Dat[0])) {
		my ($flag,$maxSC) = (0,0);
		for (@Reads) {
			my $MAPQ = $_->[4];
			my $CIGAR = $_->[5];
			while ($CIGAR =~ /(\d+)S/g) {
				$maxSC = $1 if $maxSC < $1;
			}
			$flag |= 1 if $MAPQ >= $minMAPQ;
			$flag |= 2 if $maxSC >= $minS;
		}
		if ($flag == 3) {
			for (@Reads) {
				print OUT join("\t",@$_);
			}
		}
		@Reads = ();
		$lastid = $Dat[0];
	} else {
		$lastid = $Dat[0];
		#push @Reads,\@Dat;
	}
	push @Reads,\@Dat;
}

close IN;
close OUT;
