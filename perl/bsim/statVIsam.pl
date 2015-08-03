#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

use FindBin 1.51 qw($RealBin);
use lib "$RealBin";
use simVirusInserts;

die "Usage: $0 <Host> <Virus> <bam sort-n> [more bam files]\n" if @ARGV <3;

my $Reff = shift;
my $Virf = shift;

#$Reff='hs_ref_GRCh38.p2_chr18.mfa.gz';
#$Virf='HBV.AJ507799.2.fa.gz';

my $Refh = openfile($Reff);
my $RefID = getRefChr1stID($Refh);
close $Refh;
my $Virfh = openfile($Virf);
my $VirID = getRefChr1stID($Virfh);
close $Virfh;
warn "[!]Ref:[$RefID], Virus:[$VirID].\n";

sub getInsertPos($$$) {
	my ($r1fr,$rSMS,$r12)=@_;
	my ($FiveSkip,$ThreeSkip,$FiveT,$ThreeT)=(0,0,0);
	if ($r1fr eq 'f') {
		$FiveSkip = $$rSMS[0];
		$ThreeSkip = $$rSMS[2];
	} elsif ($r1fr eq 'r') {
		$FiveSkip = $$rSMS[2];
		$ThreeSkip = $$rSMS[0];
	}
}
sub InsertPos2RealPos() {
	my ($t)=@_;
}
sub getRealPos($$$$$$) {
	my ($r1fr,$MappedChr,$rSMS,$r12,$strand,$samPos)=@_;
	getInsertPos($r1fr,$rSMS,$r12);
	return (1,2);
}

my %Stat;
for my $bamin (@ARGV) {
	print STDERR "[!]Reading [$bamin] ";
	my $bamfh = openfile($bamin);
	while (<$bamfh>) {
		my @dat1 = split /\t/;
		$_=<$bamfh>;
		my @dat2 = split /\t/;
		die "[x]Read1 & Read2 not match ! [$dat1[0]] ne [$dat2[0]]\n" if $dat1[0] ne $dat2[0];
		# sf0_Ref_2707868_2708068_2708268_Vir_-_5629_5731
		$dat1[0] =~ /^s([fr])(\d+)_Ref_(\d+)_(\d+)_(\d+)_Vir_([+-])_(\d+)_(\d+)_R__(\d+)_(\d+)$/ or die;
		my ($r1fr,$innerPos,$RefLeft,$RefMiddle,$RefRight,$VirStrand,$VirLeft,$VirRight,$InsertSize,$ReadLen) = ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10);
		#print "$dat1[0] $innerPos,$RefLeft,$RefMiddle,$RefRight,$VirStrand,$VirLeft,$VirRight\n";
		my $r1SMS = cigar2SMS($dat1[5]);
		my $r2SMS = cigar2SMS($dat2[5]);
		print "$dat1[0] [",join('.',@$r1SMS),"][",join('.',@$r1SMS),"]\n";
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
		($R1Left,$R1Right) = getRealPos($r1fr,$refR1,$r1SMS,$r12R1,$strandR1,$dat1[3]);
		warn "($R1Left,$R1Right) $r12R1 $r1fr";
	}
	close $bamfh;
	print STDERR ".\n";
}


ddx \%Stat;
__END__
./statVIsam.pl hs_ref_GRCh38.p2_chr18.mfa.gz HBV.AJ507799.2.fa.gz bam/*.sn.bam
