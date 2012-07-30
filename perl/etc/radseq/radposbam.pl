#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Data::Dump qw(ddx);
use Galaxy::IO;
use Galaxy::IO::FASTAQ qw(readfq getQvaluesFQ);
use Galaxy::SeqTools;

die "Usage: $0 <reference> <merged bam> <out>\n" if @ARGV<2;
# no more "<sorted bam files>", use `samtools merge -r <out.bam> <in1.bam> <in2.bam> [...]` first !
my $fafs=shift;
my $bamfs=shift;
my $outfs=shift;

my $ReadLen = 101-5;
my $minOKalnLen = 30;
my $minAlignLen = int($ReadLen * 0.6);
$minAlignLen = $minOKalnLen if $minAlignLen < $minOKalnLen;
die if $ReadLen < $minOKalnLen;

my $Eseq="CTGCAG";
my $EseqLen = length $Eseq;
my $EcutAt=5;
my $EseqL="CTGCA";#substr $Eseq,0,$EcutAt;
my $EseqR="TGCAG";

my $EfwdTerm5=1-$EcutAt+1;	# -3
my $EfwdTerm3=1-$EcutAt+$ReadLen;	# 91
my $ErevTerm3=0;	# 0
my $ErevTerm5=$ErevTerm3-$ReadLen+1;	# -94
# 12008 -> [12005,12099],[11914,12008]
my $Rfwd2Ec = -$EfwdTerm5;	# 3
my $Rrev2Ec = -$ErevTerm3;	# 0

my %Stat;
my $t;
open O,'>',$outfs or die "Error opening $outfs : $!\n";
$t = "# Ref: [$fafs], Enzyme: [$Eseq], Cut after $EcutAt\n# Bams: [$bamfs]\n\n";
print O $t;
print $t;

my $FHref = openfile($fafs);
my @aux = undef;
my ($name, $comment, $seq, %RefSeq, %RefLen, $ret);
my ($n, $slen) = (0, 0);
while ($ret = &readfq($FHref, \@aux)) {
	($name, $comment, $seq) = @$ret;
	++$n;
	$RefSeq{$name} = $seq;
	$RefLen{$name} = length($seq);
	$slen += $RefLen{$name};
	#warn "$name, $comment, $seq\n";
}
warn "Ref: $n seq of $slen bp.\n";
my $th = openpipe('samtools view -H',$bamfs);
my @RefOrder;
while (<$th>) {
	/^\@SQ\tSN:([^ ]+) LN:(\d+)$/ or next;
	push @RefOrder,$1;
	die "[x]Chr:$1 Len:$2 != $RefLen{$name}.\n" if $2 != $RefLen{$name};
}
close $th;
warn "Input files match.\n";

sub cigar2reflen($) {
	my $cigar = $_[0];
	my @cigar = $cigar =~ /(\d+)(\w)/g;
	my ($reflen,$maxM)=(0,0);
	my ($firstM,$lastM)=(0,0);
	$firstM = $cigar[0] if ($cigar[1] eq 'M');
	$lastM = $cigar[-2] if ($cigar[-1] eq 'M');
	while (@cigar) {
		my ($len,$op) = splice(@cigar,0,2);
		if ($op eq 'M') {
			$reflen += $len;
			$maxM = $len if $maxM < $len;
			next;
		}
		$reflen += $len if $op eq 'D';
	}
	return [$reflen,$maxM,$firstM,$lastM];
}
sub deal_cluster($$$) {
	my ($lastChr,$lastPos,$itemsA) = @_;
	my $mark = 0;
	if ( (not exists $RefSeq{$lastChr}) or $lastPos<2 or $lastPos+$EcutAt>$RefLen{$name} ) {
		++$Stat{'Cluster_err'} if $lastChr ne '';
		return [$mark,''];
	}
	++$Stat{'Cluster_cnt'};
	my $ref = substr $RefSeq{$lastChr},$lastPos-2,$EseqLen;
	my $refmatches = ($ref ^ $Eseq) =~ tr/\0//;	# http://perlmonks.org/?node_id=840593
	++$Stat{'Cluster_RefFullECcnt'}{$refmatches};
#warn ">$lastChr,$lastPos,$ref,$refmatches\n";
	my %strand;
	for (@$itemsA) {
		my ($Strand,$readSeq,$readQ) = @$_;
		my @Qvals = @{getQvaluesFQ($readQ)};
		++$strand{$Strand};
		my $mask = $readSeq ^ $EseqR;	# http://stackoverflow.com/questions/4709537/fast-way-to-find-difference-between-two-strings-of-equal-length-in-perl
		my $seqmatches = $mask =~ tr/\0//;
		while ($mask =~ /[^\0]/g) {
			my $qvar = $Qvals[$-[0]];
#warn substr($readSeq,$-[0],1), ' ', substr($EseqR,$-[0],1), ' ', $-[0], " $qvar\n";
			$mark += (10**($qvar/(-10)))/(3*$EcutAt);
		}
		$mark += $seqmatches/$EcutAt;
#warn "   $strandShift,$readSeq,$readQ [@Qvals] $mark\n";
	}
	my @strands = sort { $strand{$b} <=> $strand{$a} } keys %strand;
	my $flag;
	if (@strands>1 and $strand{$strands[1]} > 1) {
		$mark += 100000;
		$flag = 'A';
	} else {
		$flag = $strands[0];
	}
#warn "$mark <\n";
	++$Stat{'Cluster_mark'}{$mark};
	return [$mark,$flag];
}

$th = openpipe('samtools view -f64 -F1796',$bamfs);	# +0x40 -0x704
my ($lastChr,$lastPos) = ('',0);
my @items;
while (<$th>) {
	my ($id,$flag,$ref,$pos,$mapq,$cigar,$MRNM,$MPOS,$ISIZE,$seq,$qual,@opt) = split /\t/;
	#my @ta = split /\t/,$t;
	++$Stat{'Total_reads'};
	next unless grep /^XT:A:U$/,@opt;
	my @Samples = grep /^RG:Z:/,@opt;
	my $sample = '';
	if (@Samples == 1) {
		$sample = substr $Samples[0],5;
warn "$sample";
	}

	my ($strand,$thePos,$readSeq,$readQ,$radM)=('+');
	my ($reflen,$maxM,$firstM,$lastM)=@{cigar2reflen($cigar)};
	if ($flag & 16) {
		$strand = '-';
#		$radM = $lastM;
#	} else {
#		$radM = $firstM;
	}
#	if ($radM<$EcutAt or $maxM<$minOKalnLen or $reflen<$minAlignLen) {
	if ($maxM<$minOKalnLen or $reflen<$minAlignLen) {	# no need to check $radM after cutfq_withEC
		++$Stat{'CIAGR_skipped'}{$strand.$cigar};
#warn "$strand.$cigar $radM\n";
		next;
	}
	if ($flag & 16) {	# reverse
		$thePos = $pos + $reflen-1 + $Rrev2Ec;
		#$strand = '-';
		$readSeq = substr $seq,-$EcutAt;
		$readSeq = revcom($readSeq);
		$readQ = substr $qual,-$EcutAt;
		$readQ = reverse $readQ;
	} else {
		$thePos = $pos + $Rfwd2Ec;
		$readSeq = substr $seq,0,$EcutAt;
		$readQ = substr $qual,0,$EcutAt;
	}
	++$Stat{'CIAGR_used'}{$strand.$cigar};
	++$Stat{'Used_reads'};

	if ($ref eq $lastChr and $thePos == $lastPos) {
		#push @items,[$strandShift,$readSeq,$readQ];
	} else {
		my ($mark,$flag) = @{&deal_cluster($lastChr,$lastPos,\@items)};
		if ($mark) {
			print O join("\t",$lastChr,$lastPos,$flag,$mark,scalar @items),"\n";
#warn join("\t",$lastChr,$lastPos,$flag,$mark,scalar @items,@{$items[0]},@{$items[-1]}),"\n";
		}
		$lastChr = $ref;
		$lastPos = $thePos;
		@items = ();
		#push @items,[$strandShift,$readSeq,$readQ];
		#ddx \%Stat;
	}
	push @items,[$strand,$readSeq,$readQ,$sample];
}
ddx \%Stat;
close O;
