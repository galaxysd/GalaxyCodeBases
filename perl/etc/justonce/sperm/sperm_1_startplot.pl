#!/bin/env perl
use strict;
use warnings;
use IO::Unread qw(unread);
use Data::Dump qw(ddx);

die "Usage: $0 <input sam.gz/bam> <output>\n" if @ARGV < 2;
my ($inf,$outf)=@ARGV;

=pod
open CHR,'<','chr.lst' or die "[x]Cannot read chr.lst, use `samtools faidx` and append ChrID manually to get one!\n";
my %ChrGI2ID;
while (<CHR>) {
	my ($gi,$id) = split /\s+/,$_;
	$ChrGI2ID{$gi}=$id;
}
close CHR;
$ChrGI2ID{'='}='=';
ddx \%ChrGI2ID;
=cut

my $FH;
if ($inf =~ /\.sam\.gz$/i) {
	open $FH,'-|',"gzip -dc $inf" or die "Error opening $inf: $!\n";
} elsif ($inf =~ /\.bam$/i) {
	open $FH,'-|',"samtools view -h $inf" or die "Error opening $inf: $!\n";
} elsif ($inf =~ /\.sam$/i) {
	open $FH,'<',$inf or die "Error opening $inf: $!\n";
} else { die; }

my $line;
my %ChrLen;
while ($line = <$FH>) {
	unless ($line =~ /^\@/) {
		unread $FH, $line;
		last;
	}
	if ($line =~ /^\@SQ/) {
		my (undef,$gi,$len) = split /\s+/,$line;
		$gi =~ s/^SN://;
		$len =~ s/^LN://;
		$ChrLen{$gi} = $len;
	}
}
ddx \%ChrLen;

my %DatbyChrM;
#$DatbyChrM{$_}=[] for values %ChrGI2ID;
while ($line = <$FH>) {
	my ($id, $flag, $ref, $pos, $mapq, $CIAGR, $mref, $mpos, $isize, $seq, $qual, @OPT) = split /\t/,$line;
	#print "$id, $flag, Chr$ChrGI2ID{$ref}, $pos, $mapq, $CIAGR, Chr$ChrGI2ID{$mref}, $mpos, $isize\n";
	#my $commonChrID = $ChrGI2ID{$ref};
	my $commonChrID = $ref;
	my $posto1m = int($pos/1000000);
	++$DatbyChrM{$commonChrID}->[$posto1m];
}
#ddx \%DatbyChrM;

close $FH;

open OUT,'>',$outf or die "Error opening $outf: $!\n";
print OUT "# $inf\n";
for my $ccid (sort { "$a$b"=~/^\d+$/ ? $a<=>$b : $a cmp $b } keys %DatbyChrM) {
	print OUT "[$ccid]\n";
	my $ArrayRef = $DatbyChrM{$ccid};
	for my $i ( 0 .. $#$ArrayRef ) {
		my $t = $$ArrayRef[$i];
		print OUT (defined $t)?$t:0,"\t",1+$i,"\n";
	}
}

close OUT;

__END__
perl startplot.pl t.sam t.out
grep \\[ t.out

perl sperm_1_startplot.n.pl Blood-MAL.sort.bam Blood-MAL.nstat
perl sperm_1_startplot.n.pl Blood-MDA.sort.bam Blood-MDA.nstat
perl sperm_1_startplot.n.pl Sperm23-MDA.sort.bam Sperm23-MDA.nstat
perl sperm_1_startplot.n.pl Sperm24-MDA.sort.bam Sperm24-MDA.nstat
perl sperm_1_startplot.n.pl Sperm28-MDA.sort.bam Sperm28-MDA.nstat
perl sperm_1_startplot.n.pl SpermS01.sort.bam SpermS01.nstat
perl sperm_1_startplot.n.pl SpermS02.sort.bam SpermS02.nstat
perl sperm_1_startplot.n.pl SpermS03.sort.bam SpermS03.nstat

$ll -L *.bam *.nstat
-rw-r--r-- 1 gaoshengjie bc_tumor  31K Sep  2 00:08 Blood-MAL.nstat
-rw-r--r-- 1 gaoshengjie bc_tumor 131G Aug 12 04:37 Blood-MAL.sort.bam
-rw-r--r-- 1 gaoshengjie bc_tumor  31K Sep  2 02:50 Blood-MDA.nstat
-rw-r--r-- 1 gaoshengjie bc_tumor  56G Aug 12 18:30 Blood-MDA.sort.bam
-rw-r--r-- 1 gaoshengjie bc_tumor  29K Sep  2 03:24 Sperm23-MDA.nstat
-rw-r--r-- 1 gaoshengjie bc_tumor  20G Aug  8 11:52 Sperm23-MDA.sort.bam
-rw-r--r-- 1 gaoshengjie bc_tumor  29K Sep  2 03:56 Sperm24-MDA.nstat
-rw-r--r-- 1 gaoshengjie bc_tumor  20G Aug  8 12:10 Sperm24-MDA.sort.bam
-rw-r--r-- 1 gaoshengjie bc_tumor  29K Sep  2 04:30 Sperm28-MDA.nstat
-rw-r--r-- 1 gaoshengjie bc_tumor  20G Aug  8 12:02 Sperm28-MDA.sort.bam
-rw-r--r-- 1 gaoshengjie bc_tumor  29K Sep  2 05:03 SpermS01.nstat
-rw-r--r-- 1 gaoshengjie bc_tumor  22G Aug 11 13:06 SpermS01.sort.bam
-rw-r--r-- 1 gaoshengjie bc_tumor  29K Sep  2 05:36 SpermS02.nstat
-rw-r--r-- 1 gaoshengjie bc_tumor  22G Aug  8 12:37 SpermS02.sort.bam
-rw-r--r-- 1 gaoshengjie bc_tumor  29K Sep  2 06:05 SpermS03.nstat
-rw-r--r-- 1 gaoshengjie bc_tumor  19G Aug 11 15:12 SpermS03.sort.bam
