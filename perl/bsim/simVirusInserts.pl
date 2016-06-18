#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);
use Cwd 'abs_path';

use FindBin 1.51 qw($RealBin);
use lib "$RealBin";
use simVirusInserts;	# 同时输出甲基化与非甲基化的结果。
#use simVirusInsertsOld;	# 只输出100%甲基化的reads（四种碱基的）

my $RefNratioMax = 0.01;	# /Nn/
our $RefMratioMax = 0.02;	# masked as lower case in *.mfa.gz

die "Usage: $0 <Host> <Virus> <Outprefix> [ReadLen=90]\n" if @ARGV <3;

my ($Reff,$Virf,$outp,$ReadLen)=@ARGV;
$ReadLen = 90 unless defined $ReadLen;

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

open INI,'>',$outp.'.ini' or die $!;
my $Refabsf = abs_path($Reff);
my $Virabsf = abs_path($Virf);
print INI <<"CONTENT";
[RefFiles]
HostRef=$Refabsf
VirusRef=$Virabsf

[Output]
WorkDir=/share/users/huxs/work/bsvir/bsI
ProjectID=simVir

CONTENT

my %fqFiles;

my $ShortLen = int(0.5 + $ReadLen/2);
my $TinyLen = int(10.5 + $ReadLen/20);
my $LongLen = int(0.5 + $ReadLen*4/3);
my @PEins = (100,150,200,420);
@PEins = (150,220,350,500) if $ReadLen > 96;
@PEins = (60,80,120,250) if $ReadLen < 77;

my %Para = (
	PEinsertLen => $PEins[2],
	SeqReadLen => $ReadLen,
	RefNratioMax => $RefNratioMax,
	RefLen => $RefLen,
	VirLen => $VirLen,
	VirFrag => 2 * $ReadLen,
	OutPrefix => $outp . '_m13FG',
);
dosim($Refstr,$Virstr,\%Para);
$fqFiles{$Para{OutPrefix}} = $Para{PEinsertLen};

$Para{PEinsertLen} = $PEins[2];
$Para{VirFrag} = $ShortLen;
$Para{OutPrefix} = $outp . '_m2D';
dosim($Refstr,$Virstr,\%Para);
$fqFiles{$Para{OutPrefix}} = $Para{PEinsertLen};

$Para{PEinsertLen} = $PEins[2];
$Para{VirFrag} = $TinyLen;
$Para{OutPrefix} = $outp . '_m2Dty';
dosim($Refstr,$Virstr,\%Para);
$fqFiles{$Para{OutPrefix}} = $Para{PEinsertLen};

$Para{PEinsertLen} = $PEins[3];
$Para{VirFrag} = $LongLen;
$Para{OutPrefix} = $outp . '_m458AE';
dosim($Refstr,$Virstr,\%Para);
$fqFiles{$Para{OutPrefix}} = $Para{PEinsertLen};

$Para{PEinsertLen} = $PEins[2];
$Para{VirFrag} = $LongLen;
$Para{OutPrefix} = $outp . '_m9';
dosim($Refstr,$Virstr,\%Para);
$fqFiles{$Para{OutPrefix}} = $Para{PEinsertLen};

$Para{PEinsertLen} = $PEins[1];
$Para{VirFrag} = $LongLen;
$Para{OutPrefix} = $outp . '_m6';
dosim($Refstr,$Virstr,\%Para);
$fqFiles{$Para{OutPrefix}} = $Para{PEinsertLen};

$Para{PEinsertLen} = $PEins[1];
$Para{VirFrag} = $ShortLen;
$Para{OutPrefix} = $outp . '_m7C';
dosim($Refstr,$Virstr,\%Para);
$fqFiles{$Para{OutPrefix}} = $Para{PEinsertLen};

$Para{PEinsertLen} = $PEins[0];
$Para{VirFrag} = $ShortLen;
$Para{OutPrefix} = $outp . '_mB';
dosim($Refstr,$Virstr,\%Para);
$fqFiles{$Para{OutPrefix}} = $Para{PEinsertLen};

my @fps = sort keys %fqFiles;
print INI "[DataFiles]\n";
for my $i (0 .. $#fps) {
	print INI "F$fps[$i].1=",abs_path($fps[$i].'.1.fq'),"\n";
	print INI "F$fps[$i].2=",abs_path($fps[$i].'.2.fq'),"\n";
}
print INI "\n[InsertSizes]\n";
for my $i (0 .. $#fps) {
	print INI "F$fps[$i]=",$fqFiles{$fps[$i]},"\nF$fps[$i].SD=$i\n";
}
close INI;

__END__
./simVirusInserts.pl hs_ref_GRCh38.p2_chr18.mfa.gz X04615.fa.gz simout

grep \> simout_*.Ref.fa | sed 's/^simout_m//'|sed 's/.Ref.fa:>/\t/'|sed 's/Ref_/Ref:/g'|sed 's/Vir_/Vir:/'|sed 's/R_/R:/'|sed 's/_/ /g'|cat -n >simed.lst

