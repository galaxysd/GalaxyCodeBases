#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120330
=cut
use strict;
use warnings;
#use Data::Dump qw(ddx);

die "Usage: $0 <Cut sites> <sorted bam> <out prefix>\n" if @ARGV<2;
my $inec=shift;
my $insam=shift;
my $outp=shift;

my $ReadLen = 100-5;	# 95
my $minOKalnLen = 30;
my $minAlignLen = int($ReadLen * 0.75);
$minAlignLen = 33 if $minAlignLen < 33;
die if $ReadLen < 33;

my $Eseq="CTGCAG";
my $EcutAt=5;
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
	} else {die "[x]Only support .sam(.gz) & .bam [$filename]\n";}
    return $infile;
}

my $ecin = openfile($inec);
my $samin = opensam($insam);
open O,'|-',"gzip -9c >$outp.edep.gz" or die "Error opening $outp.edep.gz with gzip: $!\n";
open L,'>',$outp.'.edep.log' or die "Error opening $outp.edep.log: $!\n";
select(L);
$|=1;
print L "From [$samin],[$ecin] to [$outp.edep.gz]\n";

my (%eCut,%eDat,@ChrOrder);
while (<$ecin>) {
	next if /^(#|$)/;
	my ($chr,$pos) = split /\t/;
	#push @{$eCut{$chr}},$pos;
	++$eCut{$chr}{$pos};
	push @ChrOrder,$chr unless exists $eCut{$chr};
	$eDat{$chr.'\t'.$pos}=[0,0,0,0,0,0,0,0];	# [SumF,CountF, SumR,CountR] for average depth
	# first 4, Unique(XT:A:U); second 4, Repeat(XT:A:R).
	# XT:A:N is fragmental like "61S12M2D3M2D10M14S".
	# Mate-sw(XT:A:M) is either short or with many(like 7) mismatches(including 'N' on reads)
}
close $ecin;

sub cigar2reflen($) {
	my $cigar = $_[0];
	my @cigar = $cigar =~ /(\w)(\d+)/g;
	my ($reflen,$maxM)=(0,0);
	while (@cigar) {
		my ($op,$len) = splice(@cigar,0,2);
		if ($op eq 'M') {
			$reflen += $len;
			$maxM = $len if $maxM < $len;
			next;
		}
		$reflen += $len if $op eq 'D';
	}
	return [$reflen,$maxM];
}

my (%ChrLen);
while (<$samin>) {
	if (/^@\w\w\t\w\w:/) {
		print O $_;
		if (/^\@SQ\tSN:(\S+)\tLN:(\d+)$/) {
			if (exists $eCut{$1}) {
				$ChrLen{$1} = $2;
				print STDERR "Chr:[$1]\tLen:[$2], Cut:[",scalar @{$eCut{$1}},"]\n";
			} else {
				warn "Chr:[$1], Len:[$2] not cut.\n";
			}
		}
		next;
	}
	my @read1=split /\t/,$_,12;
	next if $read1[11] =~ /\bXT:A:[RN]\b/;
	my ($reflen,$maxM)=@{cigar2reflen($read1[5])};
	if ($maxM<$minOKalnLen or $reflen<$minAlignLen) {
		next;
	}
	my $thePos;
	if ($read1[1] & 16) {	# reverse
		$thePos = $read1[3] + $reflen-1 + $Rrev2Ec;
	} else {
		$thePos = $read1[3] + $Rfwd2Ec;
	}
}

close $samin;
close O;
print L "\n";
close L;
