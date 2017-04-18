#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);
use Cwd qw(abs_path getcwd);

use FindBin 1.51 qw($RealBin);
use lib "$RealBin";
use simVirusInserts;	# 同时输出甲基化与非甲基化的结果。
#use simVirusInsertsOld;	# 只输出100%甲基化的reads（四种碱基的）

my $RefNratioMax = 0.01;	# /Nn/
our $RefMratioMax = 0.02;	# masked as lower case in *.mfa.gz
($RefNratioMax,$RefMratioMax)=(0,0);	# 现在需要更好的结果，です！

die "Simulate Directional libraries of Bisulfite-Sequencing data.\nUsage: $0 <Host> <Virus> <Outprefix> [ReadLen=90 [Ticks.ini]]\nInvoke as: mkdir sim90 && cd sim90 && $0 && cd ..\n" if @ARGV <3;

my ($Reff,$Virf,$outp,$ReadLen,$TicksINI)=@ARGV;
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
my $Cwdabs = abs_path(getcwd()."/../run");
print INI <<"CONTENT";
[RefFiles]
HostRef=$Refabsf
VirusRef=$Virabsf

[Output]
WorkDir=$Cwdabs
ProjectID=$outp

CONTENT

my %fqFiles;

my $ShortLen = int(0.5 + $ReadLen/2);
my $TinyLen1 = 5;
my $TinyLen2 = 10;
my $TinyLen3 = 20;
my $LongLen = int(0.5 + $ReadLen*4/3);
my @PEins = (100,150,200,420);
my $maxPEins = 750;	# [77,96]
if ($ReadLen > 96) {
	@PEins = (150,220,350,530);
	$maxPEins = 500;
}
if ($ReadLen < 77) {
	@PEins = (60,80,120,250);
	$maxPEins = 250;
}
warn "ReadLen:$ReadLen, PE_ins:@PEins, maxPEins:$maxPEins.\n";

my %Para = (
	PEinsertLen => $PEins[2],
	SeqReadLen => $ReadLen,
	RefNratioMax => $RefNratioMax,
	RefLen => $RefLen,
	VirLen => $VirLen,
	VirFrag => 2 * $ReadLen,
	OutPrefix => $outp . '_m13FG',
);

my ($pRefticks,$pVirticks);
if (defined $TicksINI and -f $TicksINI) {
	print STDERR "Loading Ticks from [$TicksINI]: ";
	my $inih = openfile($TicksINI);
	while (<$inih>) {
		chomp;
		if (/^Refticks=\d/) {
			my @dat = split /=/,$_;
			@dat = split /,/,$dat[1];
			$pRefticks = \@dat;
			print STDERR "Ref:",scalar @dat,', ';
		} elsif (/^Virticks=\d/) {
			my @dat = split /=/,$_;
			@dat = split /,/,$dat[1];
			$pVirticks = \@dat;
			print STDERR "Vir:",scalar @dat,', ';
		}
	}
	my $virtickmax = $#$pVirticks;
	if ($virtickmax < $#$pRefticks) {
		for my $i ( (1+$virtickmax) .. $#$pRefticks ) {
			my $j = int($i - 0.8*$virtickmax);	# must use `$virtickmax` instead of `$#$pVirticks` !
			while ( abs($j)>$virtickmax ) {
				$j = abs($j) - $virtickmax;
			}
			$pVirticks->[$i] = $pVirticks->[$j];
		}
		print STDERR "Vir padded. ";
	}
	print STDERR "\b\b.\n";
} else {
	my $SeqReadLen = $Para{SeqReadLen};
	my $RefBorder = $maxPEins + 1000;
	$pRefticks = getticks($RefBorder,$Refstr,$RefLen,$maxPEins,$RefNratioMax);
	$pVirticks = getticks($Para{VirFrag},$Virstr,$VirLen,int(0.9+ 0.5*$Para{VirFrag}),$RefNratioMax);
}
$Para{pRefticks} = $pRefticks;
$Para{pVirticks} = $pVirticks;

dosim($Refstr,$Virstr,\%Para);
$fqFiles{$Para{OutPrefix}} = [$Para{PEinsertLen},$Para{VirFrag}];

$Para{PEinsertLen} = $PEins[2];
$Para{VirFrag} = $ShortLen;
$Para{OutPrefix} = $outp . '_m2D';
dosim($Refstr,$Virstr,\%Para);
$fqFiles{$Para{OutPrefix}} = [$Para{PEinsertLen},$Para{VirFrag}];

$Para{PEinsertLen} = $PEins[2];
$Para{VirFrag} = $TinyLen1;
$Para{OutPrefix} = $outp . '_m2Dty1';
dosim($Refstr,$Virstr,\%Para);
$fqFiles{$Para{OutPrefix}} = [$Para{PEinsertLen},$Para{VirFrag}];

$Para{PEinsertLen} = $PEins[2];
$Para{VirFrag} = $TinyLen2;
$Para{OutPrefix} = $outp . '_m2Dty2';
dosim($Refstr,$Virstr,\%Para);
$fqFiles{$Para{OutPrefix}} = [$Para{PEinsertLen},$Para{VirFrag}];

$Para{PEinsertLen} = $PEins[2];
$Para{VirFrag} = $TinyLen3;
$Para{OutPrefix} = $outp . '_m2Dty3';
dosim($Refstr,$Virstr,\%Para);
$fqFiles{$Para{OutPrefix}} = [$Para{PEinsertLen},$Para{VirFrag}];

$Para{PEinsertLen} = $PEins[3];
$Para{VirFrag} = $LongLen;
$Para{OutPrefix} = $outp . '_m458AE';
dosim($Refstr,$Virstr,\%Para);
$fqFiles{$Para{OutPrefix}} = [$Para{PEinsertLen},$Para{VirFrag}];

$Para{PEinsertLen} = $PEins[2];
$Para{VirFrag} = $LongLen;
$Para{OutPrefix} = $outp . '_m9';
dosim($Refstr,$Virstr,\%Para);
$fqFiles{$Para{OutPrefix}} = [$Para{PEinsertLen},$Para{VirFrag}];

$Para{PEinsertLen} = $PEins[1];
$Para{VirFrag} = $LongLen;
$Para{OutPrefix} = $outp . '_m6';
dosim($Refstr,$Virstr,\%Para);
$fqFiles{$Para{OutPrefix}} = [$Para{PEinsertLen},$Para{VirFrag}];

$Para{PEinsertLen} = $PEins[1];
$Para{VirFrag} = $ShortLen;
$Para{OutPrefix} = $outp . '_m7C';
dosim($Refstr,$Virstr,\%Para);
$fqFiles{$Para{OutPrefix}} = [$Para{PEinsertLen},$Para{VirFrag}];

$Para{PEinsertLen} = $PEins[0];
$Para{VirFrag} = $ShortLen;
$Para{OutPrefix} = $outp . '_mB';
dosim($Refstr,$Virstr,\%Para);
$fqFiles{$Para{OutPrefix}} = [$Para{PEinsertLen},$Para{VirFrag}];

my @fps = sort keys %fqFiles;
print INI "[DataFiles]\n";
for my $i (0 .. $#fps) {
	print INI "F$fps[$i].1=",abs_path($fps[$i].'.b1.fq.gz'),"\n";
	print INI "F$fps[$i].2=",abs_path($fps[$i].'.b2.fq.gz'),"\n";
}
print INI "[Base4Files]\n";
for my $i (0 .. $#fps) {
	print INI "F$fps[$i].1=",abs_path($fps[$i].'.r1.fq.gz'),"\n";
	print INI "F$fps[$i].2=",abs_path($fps[$i].'.r2.fq.gz'),"\n";
}
print INI "\n[InsertSizes]\n";
for my $i (0 .. $#fps) {
	print INI "F$fps[$i]=",$fqFiles{$fps[$i]}->[0],"\nF$fps[$i].SD=",$i+1,"\n";
}
print INI "\n[Simed]\n";
for my $i (0 .. $#fps) {
	print INI "F$fps[$i].VirFrag=",$fqFiles{$fps[$i]}->[1],"\n";
}
print INI 'Refticks=',join(',',@$pRefticks),"\n";
print INI 'Virticks=',join(',',@$pVirticks),"\n";
close INI;

__END__
./simVirusInserts.pl hs_ref_GRCh38.p2_chr18.mfa.gz X04615.fa.gz simout

grep -h \> *_m*.Ref.fa | sed 's/^simout_m//'|sed 's/.Ref.fa:>/\t/'|sed 's/Ref_/Ref:/g'|sed 's/Vir_/Vir:/'|sed 's/R_/R:/'|sed 's/_/ /g'|cat -n >simed.lst

~/git/toGit/perl/bsim/simVirusInserts.pl ~/nas/Ref/GRCh38_no_alt_analysis_set.fna.gz ~/nas/Ref/HBV.X04615.fa.gz s150 150

grep -h \> ../../sim90/*.Ref.fa |sed 's/>Ref_/Ref:/g'|sed 's/Vir_/Vir:/'|sed 's/R_/R:/'|sed 's/_/ /g'|cat -n
