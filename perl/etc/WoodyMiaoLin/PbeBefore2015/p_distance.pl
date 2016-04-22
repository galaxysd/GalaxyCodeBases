#!/usr/bin/perl
use strict;
use warnings;
use threads;


open IN, "<", "sequence.input";
open OUT, ">", "p_distance.log";

my %seqs;
while (<IN>) {
	next if $_ eq "\n";
	chomp;
	my @a = split /\t/;
	if ($a[0] eq "FCA6.2") {
		$seqs{FCA62} .= $a[1];
	} elsif ($a[0] eq "PBE084") {
		$seqs{PBE084} .=$a[1];
	} elsif ($a[0] eq "PBE144") {
		$seqs{PBE144} .=$a[1];
	} elsif ($a[0] eq "PVI033") {
		$seqs{PVI033} .=$a[1];
	} else {
		next;
	}
}
close IN;
warn "Read sequences complete!\n";

my %thread_compare;
$thread_compare{FCA62_PBE084} = threads -> create(\&p_distance, $seqs{FCA62}, $seqs{PBE084});
$thread_compare{FCA62_PBE144} = threads -> create(\&p_distance, $seqs{FCA62}, $seqs{PBE144});
$thread_compare{FCA62_PVI033} = threads -> create(\&p_distance, $seqs{FCA62}, $seqs{PVI033});
$thread_compare{PBE084_PBE144} = threads -> create(\&p_distance, $seqs{PBE084}, $seqs{PBE144});
$thread_compare{PBE084_PVI033} = threads -> create(\&p_distance, $seqs{PBE084}, $seqs{PVI033});
$thread_compare{PBE144_PVI033} = threads -> create(\&p_distance, $seqs{PBE144}, $seqs{PVI033});
warn "Comparing sequences...\n";

my %distance;
foreach (keys %thread_compare) {
	$distance{$_} = $thread_compare{$_} -> join();
}
warn "Compare sequences complete!\n";

my %thread_hetero;
foreach (keys %seqs) {
	$thread_hetero{$_} = threads -> create (\&heterozygous, $seqs{$_});
}

my %titv;
foreach (keys %seqs) {
	$titv{$_} = $thread_hetero{$_} -> join();
}
warn "Calculate heterozygous sites complete!\n";

print OUT "#SamplePair\t#ZeroSame\t#OneSame\t#TwoSame\t#NotN\n"; 
foreach (sort keys %distance) {
	print OUT $_;
	foreach (@{$distance{$_}}) {
		print OUT "\t$_";
	}
	print OUT "\n";
}
print OUT "#Sample\t#Ti\t#Tv\t#NotN\n";
foreach (sort keys %titv) {
	print OUT "$_\t$titv{$_}[0]\t$titv{$_}[1]\t$titv{$_}[2]\n";
}
close OUT;
warn "Done!\n";

sub p_distance {
	my ($seq1, $seq2) = @_;
	my $len1 = length $seq1;
	my $len2 = length $seq2;
	die "Different sequence length!" if $len1 != $len2;
	my $two; # both alleles are same.
	my $one; # only one allele is same.
	my $zero; # both alleles are different.
	my $total; # total informative bp;
	foreach my $i (0 .. $len1 - 1) {
		my $s1 = substr $seq1, $i, 1;
		my $s2 = substr $seq2, $i, 1;
		if (($s1 eq "N") or ($s2 eq "N")) {
			next;
		} else {
			++$total;
			if ($s1 eq $s2) {
				++$two;
			} elsif ((($s1 eq "A") and ($s2 =~ /[WMR]/))
			or (($s1 eq "C") and ($s2 =~ /[SMY]/))
			or (($s1 eq "G") and ($s2 =~ /[SKR]/))
			or (($s1 eq "T") and ($s2 =~ /[WKY]/))
			or (($s1 eq "W") and ($s2 =~ /[ATMKRY]/))
			or (($s1 eq "S") and ($s2 =~ /[CGMKRY]/))
			or (($s1 eq "M") and ($s2 =~ /[ACWSRY]/))
			or (($s1 eq "K") and ($s2 =~ /[GTWSRY]/))
			or (($s1 eq "R") and ($s2 =~ /[AGWSMK]/))
			or (($s1 eq "Y") and ($s2 =~ /[CTWSMK]/))) {
				++$one;
			} else {
				++$zero;
			}
		}
	}
	return [$zero, $one, $two, $total];
}

sub heterozygous {
	my $seq = $_[0];
	my $ti = ($seq =~ tr/[RY]//);
	my $tv = ($seq =~ tr/[WSMK]//);
	my $homo = ($seq =~ tr/AGCT//);
	my $nn = $ti + $tv + $homo;
	return [$ti, $tv, $nn];
}

