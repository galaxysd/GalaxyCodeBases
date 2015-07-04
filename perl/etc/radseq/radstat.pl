#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120330
=cut
use strict;
use warnings;
#use Data::Dump qw(ddx);

die "Usage: $0 <Cut sites> <sorted bam> <out prefix>\n" if @ARGV<3;
my $inec=shift;
my $insam=shift;
my $outp=shift;

my $ReadLen = 100-5;	# 95
my $minOKalnLen = 30;
my $minAlignLen = int($ReadLen * 0.6);
$minAlignLen = $minOKalnLen if $minAlignLen < $minOKalnLen;
die if $ReadLen < $minOKalnLen;

my $Eseq="CTGCAG";
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

sub openfile($) {
    my ($filename)=@_;
    my $infile;
    if ($filename=~/.bz2$/) {
	    open( $infile,"-|","bzip2 -dc $filename") or die "Error opening $filename: $!\n";
    } elsif ($filename=~/.gz$/) {
     	open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
    } else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
    return $infile;
}
sub opensam($) {
    my ($filename)=@_;
    my $infile;
    if ($filename=~/.bam$/) {
	    open( $infile,"-|","samtools view -h -F 132 $filename") or die "Error opening $filename: $!\n";
    } elsif ($filename=~/.sam.gz$/) {
     	open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
    } elsif ($filename=~/.sam$/) {
     	open( $infile,"<",$filename) or die "Error opening $filename: $!\n";
	} else {
		#die "[x]Only support .sam(.gz) & .bam [$filename]\n";
		open( $infile,"-|","samtools view -h -F 132 $filename") or die "Error opening $filename: $!\n. Only support .sam(.gz) & .bam !";
	}
    return $infile;
}

my $ecin = openfile($inec);
my $samin = opensam($insam);
my ($O,$N);
open $O,'|-',"gzip -9c >$outp.elen.gz" or die "Error opening $outp.edep.gz with gzip: $!\n";
open $N,'|-',"gzip -9c >$outp.nlen.gz" or die "Error opening $outp.ndep.gz with gzip: $!\n";
open L,'>',$outp.'.eclog' or die "Error opening $outp.edep.log: $!\n";
select(L);
$|=1;
select(STDOUT);
my $tmpstr="From [$insam],[$inec] to [$outp.{e,n}dep.gz]\n";
print L $tmpstr;
warn $tmpstr;
no warnings;
$tmpstr = "#".join("\t",qw/ChrID eCut_pos isMutSite Reads_Count Reads_Len_mean Reads_Cnt(fwd_U,rev_U,fwd_R,rev_R) Reads_Len_mean(same)/)."\n";
use warnings;
print $O $tmpstr;
print $N $tmpstr;

my (%nDat,%eDat,@ChrOrder);
my ($eCntAll,%eCnt)=(0);
while (<$ecin>) {
	next if /^(#|$)/;
	chomp;
	my ($chr,$pos) = split /\t/;
	unless (exists $eDat{$chr}) {
		push @ChrOrder,$chr;
		$nDat{$chr} = {};
	}
	++$eCntAll;
	++$eCnt{$chr};
	$eDat{$chr}{$pos}=[0,0,0,0,0,0,0,0,0];	# [SumF,CountF, SumR,CountR] for average read length (len on ref). The last(9th) is flag for mutated eCut_site
	# first 4, Unique(XT:A:U); second 4, Repeat(XT:A:R).
	# XT:A:N is fragmental like "61S12M2D3M2D10M14S".
	# Mate-sw(XT:A:M) is either short or with many(like 7) mismatches(including 'N' on reads)
}
close $ecin;

sub cigar2reflen($) {
	my $cigar = $_[0];
	my @cigar = $cigar =~ /(\d+)(\w)/g;
	my ($reflen,$maxM)=(0,0);
	while (@cigar) {
		my ($len,$op) = splice(@cigar,0,2);
		if ($op eq 'M') {
			$reflen += $len;
			$maxM = $len if $maxM < $len;
			next;
		}
		$reflen += $len if $op eq 'D';
	}
	return [$reflen,$maxM];
}

my ($ReadsCnt,$ReadsSkipped,%eCutSiteCnt,%ChrLen)=(0,0);
print L "[ChrNFO]\nbyChr = <<TABLE\n#ChrID\tLen\teCut_Count\n";
while (<$samin>) {
	if (/^@\w\w\t\w\w:/) {
		if (/^\@SQ\tSN:(\S+)\tLN:(\d+)$/) {
			if (exists $eDat{$1}) {
				$ChrLen{$1} = $2;
				#print STDERR "Chr:[$1]\tLen:[$2], Cut:[",scalar keys %{$eDat{$1}},"]\n";
				print L "$1\t$2\t",$eCnt{$1},"\n";
			} else {
				warn "Chr:[$1], Len:[$2] not cut.\n";
				print L "$1\t$2\t0\n";
			}
		}
		next;
	}
	my @read1=split /\t/,$_,12;
	++$ReadsCnt;
	my $AmbShift=0;
	if ($read1[11] =~ /\bXT:A:[RN]\b/) {
		$AmbShift = 4;
	}
	my ($reflen,$maxM)=@{cigar2reflen($read1[5])};
	if ($maxM<$minOKalnLen or $reflen<$minAlignLen) {
		++$ReadsSkipped;
#########
		#print join("\t",@read1[0..8]),"\n**>\t$reflen,$maxM\n";
#########
		next;
	}
	my ($strandShift,$thePos)=(0);
	if ($read1[1] & 16) {	# reverse
		$thePos = $read1[3] + $reflen-1 + $Rrev2Ec;
		$strandShift = 2;
	} else {
		$thePos = $read1[3] + $Rfwd2Ec;
	}
	if (exists $eDat{$read1[2]}{$thePos}) {
		$eDat{$read1[2]}{$thePos}->[0+$strandShift+$AmbShift] += $reflen;
		++$eDat{$read1[2]}{$thePos}->[1+$strandShift+$AmbShift];
	} else {
		$nDat{$read1[2]}{$thePos} = [0,0,0,0,0,0,0,0,0] unless (exists $nDat{$read1[2]}{$thePos});
		$nDat{$read1[2]}{$thePos}->[0+$strandShift+$AmbShift] += $reflen;
		++$nDat{$read1[2]}{$thePos}->[1+$strandShift+$AmbShift];
		if (($strandShift == 0 and $read1[5] =~ /^(\d+)M/ and $1>=$EcutAt and $read1[9] =~ /^$EseqR/) or
		    ($strandShift == 2 and $read1[5] =~ /(\d+)M$/ and $1>=$EcutAt and $read1[9] =~ /$EseqL$/))
		{
			$nDat{$read1[2]}{$thePos}->[8]=1;
			$eDat{$read1[2]}{$thePos} = $nDat{$read1[2]}{$thePos};
			delete $nDat{$read1[2]}{$thePos};
		}
	}
#########
	#print join("\t",@read1[0..8]),"\n-->$thePos\t$reflen,$maxM\t$strandShift+$AmbShift\n";
#########
}

sub printDat($$$$) {
	my ($chr,$theDatHp,$fh,$isEcutSite)=@_;
	my $FPsites=0;
	for my $pos (sort {$a<=>$b} keys %{$theDatHp}) {
		my ($sA,$sB,@tA,@tB)=(0,0);
		for my $id (1,3,5,7) {
			my $tmp = $theDatHp->{$pos}->[$id];
			push @tA,$tmp;	# n
			$sA += $tmp;
			if ($tmp) {
				my $tmpval = $theDatHp->{$pos}->[$id-1];
				push @tB,int(0.5+100*$tmpval/$tmp)/100;	# avgLen
				$sB += $tmpval;
			} else {
				push @tB,'.';
			}
		}
		if ($sA) {
			$sB = int(0.5+100*$sB/$sA)/100;
			if ($theDatHp->{$pos}->[8]) {
				if ($theDatHp->{$pos}->[1]+$theDatHp->{$pos}->[3]) {
					++$eCutSiteCnt{"TP_Mut"};
				} else {
					++$eCutSiteCnt{"TP_Mut_p"};
					$theDatHp->{$pos}->[8] = 2;
				}
				++$FPsites;	# Well, it WAS FP ...
			} else {
				++$eCutSiteCnt{"${isEcutSite}P"};
				++$FPsites if ${isEcutSite} eq 'F';
			}
		} elsif (${isEcutSite} eq 'T') {
			++$eCutSiteCnt{"FN"};
		}
		print $fh join("\t",$chr,$pos,$theDatHp->{$pos}->[8],$sA,$sB,join(",",@tA),join(",",@tB)),"\n";
	}
	$eCutSiteCnt{"TN"} -= $FPsites;
}

print L "TABLE\n\n";

for my $chr (@ChrOrder) {
	&printDat($chr,$eDat{$chr},$O,'T') if exists $eDat{$chr};
	&printDat($chr,$nDat{$chr},$N,'F') if exists $nDat{$chr};
	$eCutSiteCnt{"TN"} += $ChrLen{$chr} - $eCnt{$chr};
}

close $samin;
close $O;
close $N;
print L <<CONTENT;
[Read1_stat]
Total = $ReadsCnt
Skipped = $ReadsSkipped

[eCut_stat]
Total = $eCntAll
CONTENT
for my $t (sort {$b cmp $a} keys %eCutSiteCnt) {
	print L $t,' = ',$eCutSiteCnt{$t},"\n";
}
close L;
