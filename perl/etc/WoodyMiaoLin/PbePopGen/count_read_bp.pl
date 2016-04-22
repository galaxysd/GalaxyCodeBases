#!/usr/bin/perl
use strict;
use warnings;
use threads;

open I, "-|", "ls -1 /bak/archive/projects/LeopardCat/2.BQSR_round1_bam/*.bam";
open O, ">", "3.BQSR1/count_BQSR1_bam_read_bp.txt";

my %file;
while (<I>) {
	chomp;
	$file{$_} = openfile($_);
}
close I;
print STDERR "Read file list complete!\n";

my %thread;
foreach (sort keys %file) {
	$thread{$_} = threads->new(\&countsam, $file{$_});
}

my %count;
foreach (sort keys %thread) {
	$count{$_} = $thread{$_}->join;
	print STDERR "Count $_ complete!\n";
}

print O "#File\t#Reads_all\t#MappedBase_all\t#Reads_auto\t#MappedBase_auto\t#Reads_chrX\t#MappedBase_chrX\t#Reads_chrMT\t#MappedBase_chrMT\n";
foreach (sort keys %count) {
	print O $_, "\t", (join "\t", @{$count{$_}}), "\n";
}
close O;
print STDERR "Done!\n";

sub openfile {
	my $filename = shift;
	my $infile;
	open $infile, "-|", "samtools view $filename" or die "Error opening $filename: $!\n";
	return $infile;
}

sub countsam {
	my $in = shift;
	my ($read_all, $bp_all, $read_auto, $bp_auto, $read_chrx, $bp_chrx, $read_mt, $bp_mt);
	while (<$in>) {
		my @a = split /\t/;
		my @b = $a[5] =~ /(\d+)[M=X]/g;
		my $bp = 0;
		$bp += $_ foreach @b;
		$bp_all += $bp;
		$read_all++;
		if ($a[2] eq "chrX") {
			$bp_chrx += $bp;
			$read_chrx++;
		} elsif ($a[2] eq "chrMT") {
			$bp_mt += $bp;
			$read_mt++;
		} else {
			$bp_auto += $bp;
			$read_auto++;
		}
	}
	close $in;
	return [$read_all, $bp_all, $read_auto, $bp_auto, $read_chrx, $bp_chrx, $read_mt, $bp_mt];
}
