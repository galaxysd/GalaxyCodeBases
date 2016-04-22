#!/usr/bin/perl
use strict;
use warnings;
use threads;
use File::Basename;

my @list;
open $list[0], "-|", "ls -1 ../f1.BSNP_fasta_BQSR1/*.fa.gz";

open O, ">", "3.BQSR1/count_BSNP_fasta_BQSR1.txt";

my %length = (
	chrA1	=>	240380223,
	chrA2	=>	168638799,
	chrA3	=>	140925898,
	chrB1	=>	206538554,
	chrB2	=>	152998503,
	chrB3	=>	148068395,
	chrB4	=>	142431058,
	chrC1	=>	222198629,
	chrC2	=>	159252932,
	chrD1	=>	115468741,
	chrD2	=>	88096124,
	chrD3	=>	94101111,
	chrD4	=>	94492513,
	chrE1	=>	61081816,
	chrE2	=>	61960243,
	chrE3	=>	41224383,
	chrF1	=>	70119229,
	chrF2	=>	83953389,
	chrX	=>	127282370,
	chrMT	=>	17009,
);

my %sex = (
	FCAP0072	=>	"female",
	PBEP0005	=>	"male",
	PBEP0008	=>	"female",
	PBEP0009	=>	"male",
	PBEP0010	=>	"male",
	PBEP0013	=>	"female",
	PBEP0014	=>	"female",
	PBEP0023	=>	"male",
	PBEP0025	=>	"male",
	PBEP0026	=>	"female",
	PBEP0028	=>	"male",
	PBEP0036	=>	"male",
	PBEP0038	=>	"female",
	PBEP0039	=>	"male",
	PBEP0064	=>	"female",
	PBEP0065	=>	"male",
	PBEP0066	=>	"female",
	PBEP0067	=>	"male",
	PBEP0068	=>	"male",
	PBEP0069	=>	"male",
	PVIP0012	=>	"male",
);

print O "#File\tSex\tChromosome\t%NonN\t%GC\t%Ti\t%Tv\tTi/Tv\t%SNP\t#A\t#M\t#R\t#W\t#C\t#S\t#Y\t#G\t#K\t#T\t#N\t#ACGT\t#Ti\t#Tv\t#NonN\t#Total\n";
foreach my $l (@list) {
	my %file;
	while (<$l>) {
		chomp;
		$file{$_} = openfile($_);
	}
	close $l;
	print STDERR "Read file list complete!\n";

	my %thread;
	foreach (sort keys %file) {
		$thread{$_} = threads->new(\&countfasta, $file{$_});
	}

	my %count;
	foreach (sort keys %thread) {
		$count{$_} = $thread{$_}->join;
		print STDERR "Count $_ complete!\n";
	}

	foreach my $f (sort keys %count) {
		foreach my $c (sort keys %{$count{$f}}) {
			my $g = basename $f;
			$g =~ s/.fa.gz$//;
			$g =~ /^(\w+)/;
			print O "$g\t$sex{$1}\t$c\t", join("\t", @{$count{$f}{$c}}[16 .. 21], @{$count{$f}{$c}}[0 .. 15]), "\n";
		}
	}
}
close O;

sub openfile {
	my $filename = shift;
	my $infile;
	open $infile, "-|", "zcat $filename" or die "Error opening $filename: $!\n";
	return $infile;
}

sub countfasta {
	my $in = shift;
	my %chrcount;
	while (<$in>) {
		chomp;
		s/>//;
		my $seq = <$in>;
		chomp $seq;
		my $A = ($seq =~ tr/A//); #AA
		my $M = ($seq =~ tr/M//); #AC
		my $R = ($seq =~ tr/R//); #AG
		my $W = ($seq =~ tr/W//); #AT
		my $C = ($seq =~ tr/C//); #CC
		my $S = ($seq =~ tr/S//); #CG
		my $Y = ($seq =~ tr/Y//); #CT
		my $G = ($seq =~ tr/G//); #GG
		my $K = ($seq =~ tr/K//); #GT
		my $T = ($seq =~ tr/T//); #TT
		my $N = ($seq =~ tr/N//); #N
		my $ho = $A + $C + $G + $T; #homozygote
		my $ti = $R + $Y; #transition
		my $tv = $M + $W + $S + $K; #transversion
		my $nn = $ho + $ti + $tv; # not N
		my $total = $N + $nn;
		warn "$_\t$total" if $total != $length{$_};
		my $pnn = sprintf("%.2f", $nn / $total);
		my $gcc = sprintf("%.2f", (($M + $R + $Y + $K) / 2 + $C + $G + $S) / $nn); #GC content
		my $pti = sprintf("%.4f", $ti / $nn); # %transition
		my $ptv = sprintf("%.4f", $tv / $nn); # %transversion
		my $titv = sprintf("%.2f", $ti / $tv);
		my $psnp = sprintf("%.4f", ($ti + $tv) / $nn); # %SNP
		$chrcount{$_} = [$A, $M, $R, $W, $C, $S, $Y, $G, $K, $T, $N, $ho, $ti, $tv, $nn, $total, $pnn, $gcc, $pti, $ptv, $titv, $psnp];
	}
	close $in;
	my @auto;
	foreach my $c (keys %chrcount) {
		next if $c =~ /chrMT|chrX/;
		foreach my $n (0 .. 15) {
			$auto[$n] += $chrcount{$c}[$n];
		}
		delete $chrcount{$c};
	}
	$auto[16] = sprintf("%.2f", $auto[14] / $auto[15]);
	$auto[17] = sprintf("%.2f", (($auto[1] + $auto[2] + $auto[6] + $auto[8]) / 2 + $auto[4] + $auto[5] + $auto[7]) / $auto[14]);
	$auto[18] = sprintf("%.4f", $auto[12] / $auto[14]);
	$auto[19] = sprintf("%.4f", $auto[13] / $auto[14]);
	$auto[20] = sprintf("%.2f", $auto[12] / $auto[13]);
	$auto[21] = sprintf("%.4f", ($auto[12] + $auto[13]) / $auto[14]);;
	$chrcount{autosome} = \@auto;
	return \%chrcount;
}
