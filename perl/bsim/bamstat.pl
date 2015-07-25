#!/usr/bin/env perl
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );

use Data::Dump qw(ddx);

die "Usage: $0 <Host> <Virus> <bam sort-n>\n" if @ARGV <3;

my ($Reff,$Virf,$bamin)=@ARGV;

#$Reff='/share/users/huxs/git/toGit/perl/bsim/Chr1.fa';
$Reff='t.fa';
$Virf='/share/users/huxs/work/bsvir/HBV.AJ507799.2.fa';
$bamin='/share/users/huxs/work/bsvir/bsI/Test1_aln/T.sn.bam';

my $InsertSize = 200;
my $ReadLen = 90;

sub openfile($) {
	my ($filename)=@_;
	my $infile;
	if ($filename=~/.[bs]am(.gz)?$/) {
		open( $infile,"-|","samtools view -F 256 $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.gz$/) {
		open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.bz2$/) {
		open( $infile,"-|","bzip2 -dc $filename") or die "Error opening $filename: $!\n";
	} else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
	return $infile;
}

sub getRefChr1stID($) {
	my $GENOME = $_[0];
	while (<$GENOME>) {
		s/^>//;
		/^(\S+)/ or next;
		my $seqname = $1;
		return $seqname;
	}
}

my $Refh = openfile($Reff);
my $RefID = getRefChr1stID($Refh);
close $Refh;
my $Virfh = openfile($Virf);
my $VirID = getRefChr1stID($Virfh);
close $Virfh;
warn "[!]Ref:[$RefID], Virus:[$VirID].\n";

sub getHostPos($$$$$$$$) {
	my ($thePos,$RefLeft,$RefMiddle,$RefRight,$VirLeft,$VirRight,$VirStrand,$refR) = @_;
	my $VirLen = $VirRight - $VirLeft;
	my $LeftPannelen = $RefMiddle - $RefLeft;
	my $RightPannelen = $RefRight - $RefMiddle;
	my $retPos=0;
	if ($thePos <= $LeftPannelen) {
		$retPos = $RefLeft + $thePos +1;
	} elsif ($thePos <= $LeftPannelen+$VirLen) {
		my $p1 = $LeftPannelen - $thePos -1;
		my $p2 = $thePos - $LeftPannelen - $VirLen;
		if ($refR eq 'Host') {
			$retPos = ($p1 > $p2)?$p1:$p2;
		} elsif ($refR eq 'Virus') {
			if ($VirStrand eq '+') {
				$retPos = $VirLeft + $p1 +1;	# 病毒序列涉及正负链，需要考虑R1、R2是哪条链上的。此处$refR eq 'Virus'的代码错误。
			} elsif ($VirStrand eq '-') {
				$retPos = $VirRight + $p1 -$ReadLen +2;	# 调成加二是迁就。真正要修改的是调用前针对正负链设置$thePos。
			} else {die;}
		}
	} elsif ($thePos <= $LeftPannelen+$VirLen+$RightPannelen) {
		$retPos = $RefLeft + $thePos - $VirLen +1;
	} else {die;}
	return $retPos;
}
sub cigar2rpos($$) {
	my ($cigar,$lpos) = @_;
	my $reflen = 0;
	my @cigar = $cigar =~ /(\d+)(\w)/g;
	while (@cigar) {
		my ($len,$op) = splice(@cigar,0,2);
		$reflen += $len if $op eq 'M' or $op eq 'D';
	}
	return $lpos + $reflen -1;
}

my %Stat;
my $bamfh = openfile($bamin);
while (<$bamfh>) {
	my @dat1 = split /\t/;
	$_=<$bamfh>;
	my @dat2 = split /\t/;
	die "[x]Read1 & Read2 not match ! [$dat1[0]] ne [$dat2[0]]\n" if $dat1[0] ne $dat2[0];
	# sf0_Ref_2707868_2708068_2708268_Vir_-_5629_5731
	$dat1[0] =~ /^sf(\d+)_Ref_(\d+)_(\d+)_(\d+)_Vir_([+-])_(\d+)_(\d+)$/ or die;
	my ($innerPos,$RefLeft,$RefMiddle,$RefRight,$VirStrand,$VirLeft,$VirRight) = ($1,$2,$3,$4,$5,$6,$7);
	my $rpos1 = cigar2rpos($dat1[5],$dat1[3]);
	my $rpos2 = cigar2rpos($dat2[5],$dat2[3]);
	warn "[$dat1[0]]:\n\t$innerPos,$RefLeft,$RefMiddle,$RefRight,$VirStrand,$VirLeft,$VirRight\t@dat1[1,2,3,5] $rpos1 - @dat2[1,2,3,5] $rpos2\n";
	my ($r12R1,$r12R2)=(0,0);
	if (($dat1[1] & 0x40) and ($dat2[1] & 0x80) ) {
		$r12R1 = 1; $r12R2 = 2;
	} elsif (($dat1[1] & 0x80) and ($dat2[1] & 0x40) ) {
		$r12R1 = 2; $r12R2 = 1;
	} else {die 'X';}
	my ($strandR1,$strandR2)=('+','+');
	$strandR1 = '-' if $dat1[1] & 0x10;
	$strandR2 = '-' if $dat2[1] & 0x10;
	my ($refR1,$refR2)=('NA','NA');
	if ($dat1[2] eq $RefID) {
		$refR1 = 'Host';
	} elsif ($dat1[2] eq $VirID) {
		$refR1 = 'Virus';
	} else {$refR1 = "Other:$dat1[2]";}
	if ($dat2[2] eq $RefID) {
		$refR2 = 'Host';
	} elsif ($dat2[2] eq $VirID) {
		$refR2 = 'Virus';
	} else {$refR2 = "Other:$dat2[2]";}
	my ($R1Left,$R1Right,$R2Left,$R2Right)=(0,0,0,0);
	$R1Left = getHostPos($innerPos,$RefLeft,$RefMiddle,$RefRight,$VirLeft,$VirRight,$VirStrand,$refR1);
	$R1Right = getHostPos($innerPos+$ReadLen-1,$RefLeft,$RefMiddle,$RefRight,$VirLeft,$VirRight,$VirStrand,$refR1);
	$R2Left = getHostPos($innerPos+($InsertSize-$ReadLen),$RefLeft,$RefMiddle,$RefRight,$VirLeft,$VirRight,$VirStrand,$refR2);
	$R2Right = getHostPos($innerPos+$InsertSize-1,$RefLeft,$RefMiddle,$RefRight,$VirLeft,$VirRight,$VirStrand,$refR2);
	my $flag=0;
	if ($refR1 eq 'Host') {	# SamPairsHost
		$flag |=1;
		if (($R1Left == $dat1[3]) or ($R1Right == $rpos1)) {	# Fq1anchored
			$flag |=2;
			if ($dat1[5] eq "${ReadLen}M") {
				$flag |=4;
				$flag |=8;
			} elsif ($dat1[5] =~ /^(\d+)S/) {
				$flag |=4;
				$flag |=8 if $1 + $R1Left ==0;
			} elsif ($dat1[5] =~ /(\d+)S$/) {
				$flag |=4;
				$flag |=8 if $1 + $R1Right ==0;
			}
		}
	}
	if ($refR2 eq 'Host') {
		$flag |=16;
		if (($R2Left == $dat2[3]) or ($R1Right == $rpos2)) {	# Fq2anchored
			$flag |=32;
			if ($dat2[5] eq "${ReadLen}M") {
				$flag |=64;
				$flag |=128;
			} elsif ($dat2[5] =~ /^(\d+)S/) {
				$flag |=64;
				$flag |=128 if $1 + $R2Left ==0;
			} elsif ($dat2[5] =~ /(\d+)S$/) {
				$flag |=64;
				$flag |=128 if $1 + $R2Right ==0;
			}
		}
	}
	$flag |=4096 if $refR1 eq 'Virus';	# Virus are not supported yet.
	$flag |=8192 if $refR2 eq 'Virus';
	$flag = sprintf "%016b",$flag;
	warn "Flag:$flag\t$r12R1,$r12R2 $strandR1,$strandR2 $refR1,$refR2\t$R1Left,$R1Right,$R2Left,$R2Right\n";
	++$Stat{$flag};
=pod
$flag:
	0: Other,Other
	1: Host
=cut
}
close $bamfh;

ddx \%Stat;

__END__
./bamstat.pl /share/users/huxs/git/toGit/perl/bsim/Chr1.fa /share/users/huxs/work/bsvir/HBV.AJ507799.2.fa /share/users/huxs/work/bsvir/bsI/Test1_aln/T.sn.bam

zcat /bak/seqdata/genomes/HomoGRCh38/HomoGRCh38.fa.gz|head -n3500000 > Chr1.fa

# bamstat.pl:171: {
#   "0000000000000000" => 93,
#   "0000000000001111" => 18,
#   "0000000000010000" => 4,
#   "0000000000010001" => 253,
#   "0000000000010011" => 29,
#   "0000000000010111" => 136,
#   "0000000000011111" => 242,
#   "0000000000111111" => 169,
#   "0000000001110000" => 1,
#   "0000000001111111" => 796,
#   "0000000011110000" => 26,
#   "0000000011110001" => 13,
#   "0000000011110011" => 161,
#   "0000000011110111" => 922,
#   "0000000011111111" => 2182,
#   "0001000000000000" => 121,
#   "0001000000010000" => 720,
#   "0001000011110000" => 1150,
#   "0010000000000000" => 84,
#   "0010000000000001" => 145,
#   "0010000000000011" => 77,
#   "0010000000000111" => 98,
#   "0010000000001111" => 1584,
#   "0011000000000000" => 5585,
# }

# bamstat.pl:156: {
#   "0"    => 347,
#   "1"    => 253,
#   "3"    => 136,
#   "11"   => 242,
#   "15"   => 796,
#   "21"   => 13,
#   "23"   => 922,
#   "31"   => 2182,
#   "1027" => 29,
#   "1047" => 161,
#   "2063" => 169,
#   "4096" => 9359,
# }
