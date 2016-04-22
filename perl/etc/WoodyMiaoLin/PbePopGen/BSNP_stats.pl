#!/usr/bin/perl
use strict;
use warnings;

die "This program makes statistics of BSNP output.\nAuther: Woody
Usage: $0 <SNP.gz> <nSNP> <out.prefix>\n" if @ARGV < 3;

open SNP, "-|", "zcat $ARGV[0]";
open NSN, "<", $ARGV[1];
open O1, ">", "$ARGV[2].count.txt";
open O2, ">", "$ARGV[2].auto_cumul_ratio.txt";
open O3, ">", "$ARGV[2].chrx_cumul_ratio.txt";

my @chrx_nsnp;
my @chrx_homo;
my @chrx_hete;
my @chrx;

my @auto_nsnp;
my @auto_homo;
my @auto_hete;
my @auto;

<SNP>;
while (<SNP>) {
	my @a = split /\t/;
	my $m = 0; # max posterior probability
	foreach (9 .. 18) {
		if ($a[$_] > $m) {
			$m = $a[$_];
		}
	}
	my $q; #quality
	if ($m == 1) {
		$q = 6;
	} else {
		$q = int(-log(1-$m)/log(10));
	}
	my $d; #depth
	if ($a[30] >= 60) {
		$d = 60;
	} else {
		$d = $a[30];
	}
	if ($a[0] eq "chrMT") {
		next;
	} elsif ($a[0] eq "chrX") {
		++$chrx[$d][$q];
		if ($a[4] =~ /[ACTG]/) {
			++$chrx_homo[$d][$q];
		} else {
			++$chrx_hete[$d][$q];
		}
	} else {
		++$auto[$d][$q];
		if ($a[4] =~ /[ACTG]/) {
			++$auto_homo[$d][$q];
		} else {
			++$auto_hete[$d][$q];
		}
	}
}
close SNP;

<NSN>;
while (<NSN>) {
	my @a = split / +/;
	my @b = split //, $a[4];
	my @c = split //, $a[3];
	foreach (0 .. $a[2]-1) {
		my $depth = ord($c[$_]) - 33;
		my $d; #depth
		if ($depth >= 60) {
			$d = 60;
		} else {
			$d = $depth;
		}
		my $quality = int((ord($b[$_])-33)/5);
		my $q;
		if ($quality >= 6) {
			$q = 6;
		} else {
			$q = $quality;
		}
		if ($a[0] eq "chrMT") {
			next;
		} elsif ($a[0] eq "chrX") {
			++$chrx_homo[$d][$q];
			++$chrx_nsnp[$d][$q];
			++$chrx[$d][$q];
		} else {
			++$auto_homo[$d][$q];
			++$auto_nsnp[$d][$q];
			++$auto[$d][$q];
		}
	}
}
close NSN;

print O1 "#Autosome nSNP file\n#Depth\t#Quality\t#Basepair\n";
foreach my $d (0 .. @auto_nsnp-1) { #depth;
	$auto_nsnp[$d][0] = 0 unless $auto_nsnp[$d];
	foreach my $q (0 .. @{$auto_nsnp[$d]}-1) { #quality
		$auto_nsnp[$d][$q] = 0 unless $auto_nsnp[$d][$q];
		print O1 "$d\t$q\t$auto_nsnp[$d][$q]\n";
	}
}
print O1 "#Autosome homozygous\n#Depth\t#Quality\t#Basepair\n";
foreach my $d (0 .. @auto_homo-1) { #depth;
	$auto_homo[$d][0] = 0 unless $auto_homo[$d];
	foreach my $q (0 .. @{$auto_homo[$d]}-1) { #quality
		$auto_homo[$d][$q] = 0 unless $auto_homo[$d][$q];
		print O1 "$d\t$q\t$auto_homo[$d][$q]\n";
	}
}
print O1 "#Autosome heterozygous\n#Depth\t#Quality\t#Basepair\n";
foreach my $d (0 .. @auto_hete-1) { #depth;
	$auto_hete[$d][0] = 0 unless $auto_hete[$d];
	foreach my $q (0 .. @{$auto_hete[$d]}-1) { #quality
		$auto_hete[$d][$q] = 0 unless $auto_hete[$d][$q];
		print O1 "$d\t$q\t$auto_hete[$d][$q]\n";
	}
}
print O1 "#All autosome loci\n#Depth\t#Quality\t#Basepair\n";
foreach my $d (0 .. @auto-1) { #depth;
	$auto[$d][0] = 0 unless $auto[$d];
	foreach my $q (0 .. @{$auto[$d]}-1) { #quality
		$auto[$d][$q] = 0 unless $auto[$d][$q];
		print O1 "$d\t$q\t$auto[$d][$q]\n";
	}
}
print O1 "#ChrX nSNP file\n#Depth\t#Quality\t#Basepair\n";
foreach my $d (0 .. @chrx_nsnp-1) { #depth;
	$chrx_nsnp[$d][0] = 0 unless $chrx_nsnp[$d];
	foreach my $q (0 .. @{$chrx_nsnp[$d]}-1) { #quality
		$chrx_nsnp[$d][$q] = 0 unless $chrx_nsnp[$d][$q];
		print O1 "$d\t$q\t$chrx_nsnp[$d][$q]\n";
	}
}
print O1 "#ChrX homozygous\n#Depth\t#Quality\t#Basepair\n";
foreach my $d (0 .. @chrx_homo-1) { #depth;
	$chrx_homo[$d][0] = 0 unless $chrx_homo[$d];
	foreach my $q (0 .. @{$chrx_homo[$d]}-1) { #quality
		$chrx_homo[$d][$q] = 0 unless $chrx_homo[$d][$q];
		print O1 "$d\t$q\t$chrx_homo[$d][$q]\n";
	}
}
print O1 "#ChrX heterozygous\n#Depth\t#Quality\t#Basepair\n";
foreach my $d (0 .. @chrx_hete-1) { #depth;
	$chrx_hete[$d][0] = 0 unless $chrx_hete[$d];
	foreach my $q (0 .. @{$chrx_hete[$d]}-1) { #quality
		$chrx_hete[$d][$q] = 0 unless $chrx_hete[$d][$q];
		print O1 "$d\t$q\t$chrx_hete[$d][$q]\n";
	}
}
print O1 "#All chrX loci\n#Depth\t#Quality\t#Basepair\n";
foreach my $d (0 .. @chrx-1) { #depth;
	$chrx[$d][0] = 0 unless $chrx[$d];
	foreach my $q (0 .. @{$chrx[$d]}-1) { #quality
		$chrx[$d][$q] = 0 unless $chrx[$d][$q];
		print O1 "$d\t$q\t$chrx[$d][$q]\n";
	}
}
close O1;

foreach (\@auto_homo, \@auto_hete, \@auto) { 
	my @y = &distribution(@{$_});
	foreach (@y) {
		print O2 join("\t", @{$_}), "\n";
	}
}
close O2;

foreach (\@chrx_homo, \@chrx_hete, \@chrx) { 
	my @y = &distribution(@{$_});
	foreach (@y) {
		print O3 join("\t", @{$_}), "\n";
	}
}
close O3;

sub distribution {
	my @a = @_;
	my $t;
	my $n = 1;
	foreach (@a) {
		$n = @{$_} if @{$_} > $n;
		foreach (@{$_}) {
			$t += $_;
		}
	}
	my @b;
	my @c;
	foreach my $i (0 .. @a-1) {
		foreach my $j (0 .. $n-1) {
			if ($i == 0) {
				if ($j == 0) {
					$b[0][0] = sprintf("%.4f", $a[0][0]/$t);
					$c[0][0] = $a[0][0];
				} else {
					$a[0][$j] = 0 unless $a[0][$j];
					$b[0][$j] = sprintf("%.4f", ($a[0][$j]+$c[0][$j-1])/$t);
					$c[0][$j] = ($a[0][$j]+$c[0][$j-1]);
				}
			} else {
				if ($j == 0) {
					$a[$i][0] = 0 unless $a[$i][0];
					$b[$i][0] = sprintf("%.4f", ($a[$i][0]+$c[$i-1][0])/$t);
					$c[$i][0] = ($a[$i][0]+$c[$i-1][0]);
				} else {
					my $x = 0;
					my $y = 0;
					foreach (0 .. $i-1) {
						$x += $a[$_][$j] if $a[$_][$j];
					}
					foreach (0 .. $j-1) {
						$y += $a[$i][$_] if $a[$i][$_];
					}
					$a[$i][$j] = 0 unless $a[$i][$j];
					$b[$i][$j] = sprintf("%.4f", ($a[$i][$j]+$c[$i-1][$j-1]+$x+$y)/$t);
					$c[$i][$j] = ($a[$i][$j]+$c[$i-1][$j-1]+$x+$y);
				}
			}
		}
	}
	foreach (0 .. @b-1) {
		unshift @{$b[$_]}, $_;
	}
	unshift @b, ["#", (0 .. @{$b[0]}-2)];
	return @b;
}
