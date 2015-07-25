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

sub getHostPos($$$$$$$) {
	my ($thePos,$RefLeft,$RefMiddle,$RefRight,$VirLeft,$VirRight,$VirStrand,$refR) = @_;
	my $VirLen = $VirRight - $VirLeft;
	my $LeftPannelen = $RefMiddle - $RefLeft;
	my $RightPannelen = $RefRight - $RefMiddle;
	my $retPos=0;
	if ($thePos <= $LeftPannelen) {
		$retPos = $RefLeft + $thePos +1;
	} elsif ($thePos <= $LeftPannelen+$VirLen) {
		if ($refR eq 'Host') {
			my $p1 = $LeftPannelen - $thePos -1;
			my $p2 = $thePos - $LeftPannelen - $VirLen;
			$retPos = ($p1 > $p2)?$p1:$p2;
		} elsif ($refR eq 'Virus') {
			$retPos = 1;
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
	$R1Right = getHostPos($innerPos+89,$RefLeft,$RefMiddle,$RefRight,$VirLeft,$VirRight,$VirStrand,$refR1);
	$R2Left = getHostPos($innerPos+110,$RefLeft,$RefMiddle,$RefRight,$VirLeft,$VirRight,$VirStrand,$refR2);
	$R2Right = getHostPos($innerPos+199,$RefLeft,$RefMiddle,$RefRight,$VirLeft,$VirRight,$VirStrand,$refR2);
	warn "\t$r12R1,$r12R2 $strandR1,$strandR2 $refR1,$refR2\t$R1Left,$R1Right,$R2Left,$R2Right\n";
	++$Stat{'SamPairs'};
}
close $bamfh;

__END__
./bamstat.pl /share/users/huxs/git/toGit/perl/bsim/Chr1.fa /share/users/huxs/work/bsvir/HBV.AJ507799.2.fa /share/users/huxs/work/bsvir/bsI/Test1_aln/T.sn.bam

zcat /bak/seqdata/genomes/HomoGRCh38/HomoGRCh38.fa.gz|head -n3500000 > Chr1.fa
