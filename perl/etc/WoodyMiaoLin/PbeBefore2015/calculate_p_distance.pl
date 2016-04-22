#!/usr/bin/perl
use strict;
use warnings;
use threads;

open I, "<", "input_fa.lst";
open O, ">", "RRHS_p_distance.txt";

my @thread;
while (<I>) {
	chomp;
	my $a = threads -> create(\&distance, $_);
	push @thread, $a;
}
close I;
warn "Read file list complete!\n";

print O "#fasta\tFCA62_PBE084\tFCA62_PBE144\tFCA62_PVI033\tPBE084_PBE144\tPBE084_PVI033\tPBE144_PVI033\n";
foreach (@thread) {
	my $a = $_ -> join();
	print O join("\t", @{$a}), "\n";
}
close O;
warn "Done!\n";

sub distance {
	$_[0] =~ /(sequence_\d+)\.fa$/;
	my @compare;
	$compare[0] = $1;
	open my $in, "<", $_[0];
	my %seqs;
	while (<$in>) {
		chomp;
		s/>//;
		s/\.//;
		$/ = ">";
		$seqs{$_} = <$in>;
		chomp $seqs{$_};
		$seqs{$_} =~ s/\s//g;
		$/ = "\n";
	}
	close $in;
	warn "read $compare[0] complete!\n";
	$compare[1] = &p_distance($seqs{FCA62}, $seqs{PBE084});
	$compare[2] = &p_distance($seqs{FCA62}, $seqs{PBE144});
	$compare[3] = &p_distance($seqs{FCA62}, $seqs{PVI033});
	$compare[4] = &p_distance($seqs{PBE084}, $seqs{PBE144});
	$compare[5] = &p_distance($seqs{PBE084}, $seqs{PVI033});
	$compare[6] = &p_distance($seqs{PBE144}, $seqs{PVI033});
	warn "compare $compare[0] complete!\n";
	return \@compare;
}

sub p_distance {
	my ($seq1, $seq2) = @_;
	my $len1 = length $seq1;
	my $len2 = length $seq2;
	die "Different sequence length!" if $len1 != $len2;
	my $zero; # both alleles are different.
	my $total; # total informative bp;
	foreach my $i (0 .. $len1 - 1) {
		my $s1 = substr $seq1, $i, 1;
		my $s2 = substr $seq2, $i, 1;
		if (($s1 eq "N") or ($s2 eq "N")) {
			next;
		} else {
			++$total;
			if ($s1 ne $s2) {
				++$zero;
			}
		}
	}
	my $p = $zero/$total;
	return $p;
}
