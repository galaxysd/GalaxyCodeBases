#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);
use Cwd 'abs_path';

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

open INI,'>',$outp.'.ini' or die $!;
print INI <<CONTENT;
[RefFiles]
HostRef=/bak/seqdata/genomes/HomoGRCh38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
VirusRef=/share/users/huxs/work/bsvir/X04615.fa.gz

[Output]
WorkDir=/share/users/huxs/work/bsvir/bsI
ProjectID=simVir

CONTENT

my %fqFiles;

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
$fqFiles{$Para{OutPrefix}} = $Para{PEinsertLen};

$Para{PEinsertLen} = 200;
$Para{VirFrag} = 45;
$Para{OutPrefix} = $outp . '_m2D';
dosim($Refstr,$Virstr,\%Para);
$fqFiles{$Para{OutPrefix}} = $Para{PEinsertLen};

$Para{PEinsertLen} = 420;
$Para{VirFrag} = 120;
$Para{OutPrefix} = $outp . '_m458AE';
dosim($Refstr,$Virstr,\%Para);
$fqFiles{$Para{OutPrefix}} = $Para{PEinsertLen};

$Para{PEinsertLen} = 200;
$Para{VirFrag} = 120;
$Para{OutPrefix} = $outp . '_m9';
dosim($Refstr,$Virstr,\%Para);
$fqFiles{$Para{OutPrefix}} = $Para{PEinsertLen};

$Para{PEinsertLen} = 150;
$Para{VirFrag} = 120;
$Para{OutPrefix} = $outp . '_m6';
dosim($Refstr,$Virstr,\%Para);
$fqFiles{$Para{OutPrefix}} = $Para{PEinsertLen};

$Para{PEinsertLen} = 150;
$Para{VirFrag} = 45;
$Para{OutPrefix} = $outp . '_m7C';
dosim($Refstr,$Virstr,\%Para);
$fqFiles{$Para{OutPrefix}} = $Para{PEinsertLen};

$Para{PEinsertLen} = 100;
$Para{VirFrag} = 45;
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
	print INI "F$fps[$i]=",$fqFiles{$fps[$i]},"\nF$fps[$i].SD=1\n";
}
close INI;

__END__
./simVirusInserts.pl hs_ref_GRCh38.p2_chr18.mfa.gz X04615.fa.gz simout

grep \> simout_*.Ref.fa | sed 's/^simout_m//'|sed 's/.Ref.fa:>/\t/'|sed 's/Ref_/Ref:/g'|sed 's/Vir_/Vir:/'|sed 's/R_/R:/'|sed 's/_/ /g'|cat -n >simed.lst

