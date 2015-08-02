#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

use FindBin 1.51 qw($RealBin);
use lib "$RealBin";
use simVirusInserts;

my $SampleCnt = 100;
#my $Depth = 50;
my $PEinsertLen=200;
my $SeqReadLen=90;
my $RefBorder = $PEinsertLen + 1000;
my $RefNratioMax = 0.02;

die "Usage: $0 <Host> <Virus> <Outprefix>\n" if @ARGV <3;

my ($Reff,$Virf,$outp)=@ARGV;

#$Reff='hs_ref_GRCh38.p2_chr18.mfa.gz';
#$Virf='HBV.AJ507799.2.fa.gz';

my $Refh = openfile($Reff);
my $Refstr = getRefChr1st($Refh);
my $RefLen = length $Refstr;
close $Refh;
my $Virfh = openfile($Virf);
my $Virstr = getRefChr1st($Virfh);
my $VirLen = length $Virstr;
close $Virfh;
$Virstr .= $Virstr;	# circle

my %Para = (
	PEinsertLen => $PEinsertLen,
	SeqReadLen => $SeqReadLen,
	RefNratioMax => $RefNratioMax,
	RefBorder => $RefBorder,
	LeftStart => 1,
	LeftEnd => $SeqReadLen-1,
	RefLen => $RefLen,
	VirLen => $VirLen,
	VirFrag => 2*$PEinsertLen,
	OutPrefix => $outp . '_m01',
);
dosim($Refstr,$Virstr,\%Para);

$Para{VirFrag} = int($SeqReadLen*0.5);
$Para{LeftStart} = $Para{VirFrag}+1;
$Para{LeftEnd} = $SeqReadLen -1;
$Para{OutPrefix} = $outp . '_m02';
dosim($Refstr,$Virstr,\%Para);

#my $VirFragMax = 500;
#my $VirFragMin = 20;

__END__
./simVirusInserts.pl hs_ref_GRCh38.p2_chr18.mfa.gz HBV.AJ507799.2.fa.gz simout

