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

my (%fqFiles,%Para);

my @ReadLen = qw(50 90 150);
my @VirFragLens = qw(5 10 15 20 30 50 75 90 100 120 150 180 200 250 280 300 350 400 450 500 550 600);
my @PEins = qw(150 240 360 500 750);

my $maxPEins = 750;

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
	print STDERR "\b\b.\n";
} else {
	my $SeqReadLen = $ReadLen[-1];
	my $RefBorder = $maxPEins + 1000;
	$pRefticks = getticks($RefBorder,$Refstr,$RefLen,$maxPEins,$RefNratioMax);
	$pVirticks = getticks($VirFragLens[-1],$Virstr,$VirLen,int(0.9+ 0.5*$VirFragLens[-1]),$RefNratioMax);
}
#$Para{pRefticks} = $pRefticks;
#$Para{pVirticks} = $pVirticks;

for my $vf (@VirFragLens) {
	for my $rl (@ReadLen) {
		for my $pi (@PEins) {
			next if $pi < $rl *3/2;
			%Para = (
				PEinsertLen => $pi,
				SeqReadLen => $rl,
				RefNratioMax => $RefNratioMax,
				RefLen => $RefLen,
				VirLen => $VirLen,
				VirFrag => $vf,
				OutPrefix => "${outp}_i${pi}v${vf}r${rl}",
				pRefticks => $pRefticks,
				pVirticks => $pVirticks,
			);
			dosim($Refstr,$Virstr,\%Para);
			$fqFiles{$Para{OutPrefix}} = [$Para{PEinsertLen},$Para{VirFrag}];
		}
	}
}

my @fps = sort keys %fqFiles;
print INI "[DataFiles]\n";
for my $i (0 .. $#fps) {
	print INI "F$fps[$i].1=",abs_path($fps[$i].'.1.fq'),"\n";
	print INI "F$fps[$i].2=",abs_path($fps[$i].'.2.fq'),"\n";
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
