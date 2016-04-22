#!/usr/bin/perl
use strict;
use warnings;

open my $fca62, "-|", "gzip -dc Felis_catus-6.2_masked.fa.gz";
open my $pbe084, "-|", "gzip -dc PBE084_masked.fa.gz";
open my $pbe144, "-|", "gzip -dc PBE144_masked.fa.gz";
open my $pvi033, "-|", "gzip -dc PVI033_masked.fa.gz";
open OUT, ">", "exons.nex";

my %file = (
	FCA62	=>	$fca62,
	PBE084	=>	$pbe084,
	PBE144	=>	$pbe144,
	PVI033	=>	$pvi033,
);

my %chr = (
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

my %seq;
foreach (sort keys %file) {
	my $u = $_;
	$seq{$u} = &readseq($file{$u});
	warn "Read $u sequence complete!\n";
}

my %loci;
my $ncha;
foreach my $u (sort keys %chr) {
	for (my $i=0; $i < $chr{$u}; ++$i) {
		my $w = substr($seq{FCA62}{$u}, $i, 100);
		next if (length $w) < 100;
		next if $w =~ /^N/;
		my @subseq;
		foreach my $v (sort keys %seq) {
			push @subseq, substr($seq{$v}{$u}, $i, 100);
		}
		if (grep(/N{5}/, @subseq) == 0) { # find 100bp windows without /N{5}/
			push @{$loci{$u}}, \@subseq;
			$i += 99;
		}
	}
	my $w = @{$loci{$u}};
	$ncha += $w*100;
	warn "$u exons\t$ncha bp\n";
}
warn "Pick loci complete!\n";

print OUT "#NEXUS
BEGIN DATA;
Dimensions ntax=4 nchar=$ncha;
Format datatype=DNA interleave=yes;
Matrix\n";

foreach my $u (sort keys %loci) {
	my $v;
	foreach my $w (@{$loci{$u}}) {
		++$v;	
#		print OUT "${u}_$v\t4\t1000\n";
		print OUT "FCA6.2\t${$w}[0]\n";
		print OUT "PBE084\t${$w}[1]\n";
		print OUT "PBE144\t${$w}[2]\n";
		print OUT "PVI033\t${$w}[3]\n\n";
	}
}

print OUT ";\nEND;\n";
close OUT;
warn "Done!\n";

sub readseq {
	my $in = $_[0];
	my %fa;
	while (<$in>) {
		last unless /(chr\w+)/;
		my $name = $1;
		$/ = ">";
		my $sequ = <$in>;
		chomp $sequ;
		$sequ =~ s/\s//g;
		$/ = "\n";
		$fa{$name} = $sequ;
	}
	close $in;
	return \%fa;
}
