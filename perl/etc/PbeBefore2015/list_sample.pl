#!/usr/bin/perl
use strict;
use warnings;

my %sample;

open I1, "<", "mtHaplotype_sampleID.txt";
while (<I1>) {
	next unless /^Pbe/;
	chomp;
	my @a = split /\t/;
	my @b = split / /, $a[1];
	foreach (@b) {
#		if (/^K/) {
#			$_ = "Pbe".$_;
#		}
		$sample{$_}[1] = $a[0];
	}
}
close I1;

open I2, "<", "xHaplotype_sampleID.txt";
while (<I2>) {
	next unless /^Pbe/;
	chomp;
	my @a = split /\t/;
	my @b = split / /, $a[4];
	$sample{$_}[2] = $a[0] foreach @b;
}
close I2;

open I3, "<", "yHaplotype_sampleID.txt";
while (<I3>) {
	next unless /^Pbe/;
	chomp;
	my @a = split /\t/;
	my @b = split / /, $a[2];
	$a[0] =~ s/-\d+//;
	$sample{$_}[3] = $a[0] foreach @b;
}
close I3;

open I4, "<", "samples_with_autosome_genotype.txt";
my $b = <I4>;
chomp $b;
my @c = split / /, $b;
foreach (@c) {
		$sample{$_}[4] = "yes";
}
close I4;

open I0, "<", "sampleID_locale.txt";
while (<I0>) {
	chomp;
	my @a = split /\t/;
	$sample{$a[0]}[0] = $a[1] if $sample{$a[0]};
}
close I0;

open O, ">", "Pbe_sample_list.txt";
print O "#SampleID\tLocale\tMT\tX\tY\tAutosome\n";
foreach (sort keys %sample) {
	$sample{$_}[0] = "n/a" unless $sample{$_}[0];
	$sample{$_}[1] = "n/a" unless $sample{$_}[1];
	$sample{$_}[2] = "n/a" unless $sample{$_}[2];
	$sample{$_}[3] = "n/a" unless $sample{$_}[3];
	$sample{$_}[4] = "n/a" unless $sample{$_}[4];
	print O "$_\t", join("\t", @{$sample{$_}}), "\n";
}
close O;
