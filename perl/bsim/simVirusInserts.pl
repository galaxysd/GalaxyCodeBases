#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

use FindBin 1.51 qw($RealBin);
use lib "$RealBin";
use simVirusInserts;

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
	PEinsertLen => 200,
	SeqReadLen => 90,
	RefNratioMax => $RefNratioMax,
	RefLen => $RefLen,
	VirLen => $VirLen,
	VirFrag => 2 * 90,
	OutPrefix => $outp . '_m13FG',
);
dosim($Refstr,$Virstr,\%Para);

$Para{PEinsertLen} = 200;
$Para{VirFrag} = 45;
$Para{OutPrefix} = $outp . '_m2D';
dosim($Refstr,$Virstr,\%Para);

$Para{PEinsertLen} = 420;
$Para{VirFrag} = 120;
$Para{OutPrefix} = $outp . '_m458AE';
dosim($Refstr,$Virstr,\%Para);

$Para{PEinsertLen} = 200;
$Para{VirFrag} = 120;
$Para{OutPrefix} = $outp . '_m9';
dosim($Refstr,$Virstr,\%Para);

$Para{PEinsertLen} = 150;
$Para{VirFrag} = 120;
$Para{OutPrefix} = $outp . '_m6';
dosim($Refstr,$Virstr,\%Para);

$Para{PEinsertLen} = 150;
$Para{VirFrag} = 45;
$Para{OutPrefix} = $outp . '_m7C';
dosim($Refstr,$Virstr,\%Para);

$Para{PEinsertLen} = 100;
$Para{VirFrag} = 45;
$Para{OutPrefix} = $outp . '_mB';
dosim($Refstr,$Virstr,\%Para);

__END__
./simVirusInserts.pl hs_ref_GRCh38.p2_chr18.mfa.gz HBV.AJ507799.2.fa.gz simout

