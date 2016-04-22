#!/usr/bin/perl
use strict;
use warnings;
use threads;

open I, "-|", "ls -1 ../c.sort_rmdup_bam/*.bam";
open O, ">", "count_FLAG_CIGAR_result.txt";

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

print O "#File\t#FLAG\t#CIGAR\n";
foreach (sort keys %count) {
	print O $_, "\t", (join "\t", @{$count{$_}}), "\n";
}
close O;
print STDERR "Done!\n";

sub openfile {
	my $filename = shift;
	my $infile;
	open $infile, "-|", "samtools-0.1.7 view $filename" or die "Error opening $filename: $!\n";
	return $infile;
}

sub countsam {
	my $in = shift;
	my %FLAG;
	my %CIGAR;
	while (<$in>) {
		my @a = split /\t/;
		++$FLAG{$a[1]};
		my @b = $a[5] =~ /(\d+\D)/g;
		foreach (@b) {
			/(\d+)(\D)/;
			$CIGAR{$2} += $1;
		}
	}
	close $in;
	my $flag;
	my $cigar;
	$flag .= "#$_=$FLAG{$_}," foreach sort keys %FLAG;
	$cigar .= "$CIGAR{$_}$_," foreach sort keys %CIGAR;
	return [$flag, $cigar];
}
