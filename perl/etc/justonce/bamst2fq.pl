#!/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

our $MINLENGTH = 5;

die "Usage: $0 <input sam.gz/bam> <output(gzipped)>\n" if @ARGV < 2;
my ($inf,$outf)=@ARGV;

my $FH;
if ($inf =~ /\.sam\.gz$/i) {
	open $FH,'-|',"gzip -dc $inf" or die "Error opening $inf: $!\n";
} elsif ($inf =~ /\.bam$/i) {
	open $FH,'-|',"samtools view -h $inf" or die "Error opening $inf: $!\n";
} elsif ($inf =~ /\.sam$/i) {
	open $FH,'<',$inf or die "Error opening $inf: $!\n";
} else { die; }

my ($line,$t);

open OUT,'|-', "gzip -9c >$outf" or die "Error opening $outf: $!\n";

while ($line = <$FH>) {
	next if $line =~ /^\@/;
	my ($id, $flag, $ref, $pos, $mapq, $CIAGR, $mref, $mpos, $isize, $seq, $qual, @OPT) = split /\t/,$line;
	#print "$id, $flag, $ref, $pos, $mapq, $CIAGR, $mref, $mpos, $isize\n";
	$t = 0;
	my @Ss=();
	while ($CIAGR =~ /(\d+)(\D)/g) {
		if ($2 eq 'S') {
			push @Ss,[$t,$1] if $1 >= $MINLENGTH;
		}
		$t += $1 if $2 ne 'D';
	}
	$t = 0;
	for (@Ss) {
		++$t;
		#if ($_->[0] + $_->[1] -1 > length $seq) {
		#	ddx \@Ss;
		#	die "$line";
		#}
		my $outseq = substr $seq,$_->[0],$_->[1];
		my $outqual = substr $qual,$_->[0],$_->[1];
		print OUT join("\n",'@'."$id.$t $flag $CIAGR",$outseq,'+',$outqual),"\n";
	}
}

close $FH;
close OUT;
__END__

perl -e 'use Data::Dump qw(ddx); $a="3S22M7S1D9S";$t=0;@Ss=();
 while ($a=~/(\d+)(\D)/g) {print pos $a,"\t$1 $2\n";
 {push @Ss,[$t,$1,$t+$1-1,$2];} $t += $1 if $2 ne 'D';} ddx \@Ss;'

find xtubam/*bam | while read a;do echo perl bamst2fq.pl $a $a.S.fq.gz \&;done
perl bamst2fq.pl xtubam/mdaSperm23xtu.bam xtubam/mdaSperm23xtu.bam.S.fq.gz &
