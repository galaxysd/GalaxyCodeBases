#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);
use Cwd qw(abs_path getcwd);
use Parallel::ForkManager;

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

my (%fqFiles,%Para,%Merge);

my @ReadLen = qw(50 90 150);
my @VirFragLens = qw(5 10 15 20 30 50 75 90 100 120 150 180 200 250 280 300 350 400 450 500 550 600);
my @PEins = qw(250 500 750);

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

my $pm = Parallel::ForkManager->new(8);	# https://metacpan.org/pod/Parallel::ForkManager#Data-structure-retrieval
#@ReadLen = (50, 90); @VirFragLens = (10, 20); @PEins = (250, 500);

# data structure retrieval and handling
$pm -> run_on_finish ( # called BEFORE the first call to start()
	sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
		# retrieve data structure from child
		if (defined($data_structure_reference)) {  # children are not forced to send anything
			my @Dat = @{$data_structure_reference};  # child passed a array reference
			#print ":@Dat\n";
			$fqFiles{$Dat[0]} = [$Dat[1],$Dat[2]];
			push @{$Merge{$Dat[3]}{$Dat[4]}},$Dat[0];
		} else {  # problems occurring during storage or retrieval will throw a warning
			print qq|[x]No message received from child process $pid!\n|;
		}
	}
);

for my $vf (@VirFragLens) {
	for my $rl (@ReadLen) {
		DOIT:
		for my $pi (@PEins) {
			next if $pi < $rl *4/3;
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
			$pm->start and next DOIT;
			dosim($Refstr,$Virstr,\%Para);
			#$fqFiles{$Para{OutPrefix}} = [$Para{PEinsertLen},$Para{VirFrag}];
			#push @{$Merge{$pi}{$vf}},$Para{OutPrefix};
			$pm->finish(0, [$Para{OutPrefix},$Para{PEinsertLen},$Para{VirFrag},$pi,$vf]);
		}
	}
}
$pm->wait_all_children;

my @fps = sort keys %fqFiles;
print INI "[DataFiles]\n";
for my $i (0 .. $#fps) {
	print INI "F$fps[$i].1=",abs_path($fps[$i].'.1.fq.gz'),"\n";
	print INI "F$fps[$i].2=",abs_path($fps[$i].'.2.fq.gz'),"\n";
}
for my $pi (@PEins) {
	for my $vf (@VirFragLens) {
		next unless exists $Merge{$pi}{$vf};
		my @OutPrefix = @{$Merge{$pi}{$vf}};
		print INI "Mv${vf}i${pi}.1=",join(',',(map {abs_path($_.'.1.fq.gz')} @OutPrefix)),"\n";
		print INI "Mv${vf}i${pi}.2=",join(',',(map {abs_path($_.'.2.fq.gz')} @OutPrefix)),"\n";
	}
}
print INI "\n[InsertSizes]\n";
for my $i (0 .. $#fps) {
	print INI "F$fps[$i]=",$fqFiles{$fps[$i]}->[0],"\nF$fps[$i].SD=",10+$i/100,"\n";
}
my $i=0;
for my $pi (@PEins) {
	for my $vf (@VirFragLens) {
		next unless exists $Merge{$pi}{$vf};
		++$i;
		print INI "Mv${vf}i${pi}=",500,"\nMv${vf}i${pi}.SD=",$i/10,"\n";
	}
}
print INI "\n[Simed]\n";
for my $i (0 .. $#fps) {
	print INI "F$fps[$i].VirFrag=",$fqFiles{$fps[$i]}->[1],"\n";
}
for my $pi (@PEins) {
	for my $vf (@VirFragLens) {
		next unless exists $Merge{$pi}{$vf};
		print INI "Mv${vf}i${pi}.VirFrag=$vf\n";
	}
}
print INI 'Refticks=',join(',',@$pRefticks),"\n";
print INI 'Virticks=',join(',',@$pVirticks),"\n";
close INI;

__END__
