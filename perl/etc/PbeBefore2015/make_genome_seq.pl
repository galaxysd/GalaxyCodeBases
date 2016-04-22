#!/usr/bin/perl

# This program make genome sequence using reference genome and BSNP output.
use strict;
use warnings;

die "Usage: $0 <reference.fa.gz> <SNP.gz> <nSNP.gz> <output.fa.gz>\n" if @ARGV < 4;

open REF, "-|", "gzip -dc $ARGV[0]";
open SNP, "-|", "gzip -dc $ARGV[1]";
open NSN, "-|", "gzip -dc $ARGV[2]";
open OUT, "|-", "gzip -9c >$ARGV[3]";

my %chr = (
	chrA1	=>	"gi|362110686|gb|CM001378.1|",
	chrA2	=>	"gi|362110683|gb|CM001379.1|",
	chrA3	=>	"gi|362110681|gb|CM001380.1|",
	chrB1	=>	"gi|362110674|gb|CM001381.1|",
	chrB2	=>	"gi|362110670|gb|CM001382.1|",
	chrB3	=>	"gi|362110654|gb|CM001383.1|",
	chrB4	=>	"gi|362110649|gb|CM001384.1|",
	chrC1	=>	"gi|362110644|gb|CM001385.1|",
	chrC2	=>	"gi|362110642|gb|CM001386.1|",
	chrD1	=>	"gi|362110640|gb|CM001387.1|",
	chrD2	=>	"gi|362110638|gb|CM001388.1|",
	chrD3	=>	"gi|362110636|gb|CM001389.1|",
	chrD4	=>	"gi|362110634|gb|CM001390.1|",
	chrE1	=>	"gi|362110632|gb|CM001391.1|",
	chrE2	=>	"gi|362110629|gb|CM001392.1|",
	chrE3	=>	"gi|362110625|gb|CM001393.1|",
	chrF1	=>	"gi|362110622|gb|CM001394.1|",
	chrF2	=>	"gi|362110619|gb|CM001395.1|",
	chrX	=>	"gi|362110616|gb|CM001396.1|",
);
my %gi = reverse %chr;

my %chr_len = (
	chrA1	=>	239302903,
	chrA2	=>	169043629,
	chrA3	=>	142459683,
	chrB1	=>	205241052,
	chrB2	=>	154261789,
	chrB3	=>	148491654,
	chrB4	=>	144259557,
	chrC1	=>	221441202,
	chrC2	=>	157659299,
	chrD1	=>	116869131,
	chrD2	=>	89822065,
	chrD3	=>	95741729,
	chrD4	=>	96020406,
	chrE1	=>	63002102,
	chrE2	=>	64039838,
	chrE3	=>	43024555,
	chrF1	=>	68669167,
	chrF2	=>	82763536,
	chrX	=>	126427096,
);

my %out = (
	chrA1	=>	"N" x 239302903,
	chrA2	=>	"N" x 169043629,
	chrA3	=>	"N" x 142459683,
	chrB1	=>	"N" x 205241052,
	chrB2	=>	"N" x 154261789,
	chrB3	=>	"N" x 148491654,
	chrB4	=>	"N" x 144259557,
	chrC1	=>	"N" x 221441202,
	chrC2	=>	"N" x 157659299,
	chrD1	=>	"N" x 116869131,
	chrD2	=>	"N" x 89822065,
	chrD3	=>	"N" x 95741729,
	chrD4	=>	"N" x 96020406,
	chrE1	=>	"N" x 63002102,
	chrE2	=>	"N" x 64039838,
	chrE3	=>	"N" x 43024555,
	chrF1	=>	"N" x 68669167,
	chrF2	=>	"N" x 82763536,
	chrX	=>	"N" x 126427096,
);

my %ref;
while (<REF>) {
	/(chr\w+)/;
	my $name = $1;
	$/ = ">";
	my $seq = <REF>;
	chomp $seq;
	$seq =~ s/\s//g;
	$/ = "\n";
	$ref{$name} = $seq;
}
close REF;
warn "Read reference complete!\n";

while (<SNP>) {
	next if /^#/;
	my @a = split /\t/;
	next if substr($ref{$gi{$a[0]}}, $a[1], 1) eq "N";
	next if $a[-4] < 5; # Next if NumGoodReads < 5
	my @rq = split //, $a[-2]; # Read Quality Scores
	my @aq = split //, $a[-1]; # Align Quality Scores
	my $nr = @rq; # Number of Reads
	my $ec; # Effective Coverage
	foreach (0 .. $nr-1) {
		my $rq = ord($rq[$_]) - 33;
		my $aq = ord($aq[$_]) - 33;
		my $c = (1 - (10**(-$rq/10))) * (1 - (10**(-$aq/10)));
		$ec += $c;
	}
	next if $ec < 5; # Effective Coverage < 5
	foreach (9 .. 18) {
		if ($a[$_] > 0.999) { # P(genotype|data) > 0.999
			substr($out{$gi{$a[0]}}, $a[1], 1) = $a[4];
			last;
		}
	}
}
close SNP;
warn "Read SNP complete!\n";

while (<NSN>) {
	next if /^#/;
	my @a = split / +/;
	my @b = split //, $a[4];
	my @c = split //, $a[3];
	foreach (0 .. $a[2]-1) {
		next if ord($c[$_]) < 38; # Coverage < 5;
		if (ord($b[$_]) >= 48) { # P(is a SNP) <= 0.001
			my $coor = $a[1] + $_; # coordinate in SNP file is O based
			substr($out{$gi{$a[0]}}, $coor, 1) = substr($ref{$gi{$a[0]}}, $coor, 1);
		}
	}
}
close NSN;
warn "Read nSNP complete!\n";

warn "#Chr\tRef Length\tNo. filtered\tNo. Remained\n";
foreach my $u (sort keys %out) {
	print OUT ">$u\n";
	my $i;
	my $i0;
	my $i1;
	foreach my $v (0 .. $chr_len{$u}-1) {
		print OUT substr($out{$u}, $v, 1);
		++$i;
		if (substr($out{$u}, $v, 1) eq "N") {
			++$i0;
		} else {
			++$i1;
		}
		if ($i == 100) {
			print OUT "\n";
			$i = 0;
		}
	}
	print OUT "\n";
	warn "$u\t$chr_len{$u}\t$i0\t$i1\n";
}
close OUT;
warn "Done!\n";
