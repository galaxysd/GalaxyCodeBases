#!/usr/bin/perl
use strict;
use warnings;

die("
Program: IndivID
Version: 1.0
Release: Mar. 6, 2013\n
Auther: Woody
Consultant: Galaxy\n
Usage: $0 <input> <output>\n
The input file should look like:\n
#ID	Locus1		Locus2		Locus3		...
Sample1	160	162	166	166	0	0	...
Sample2	0	0	174	176	178	178	...
Sample3	180	182	186	186	188	188	...
...	...	...	...	...	...	...	...
\nMissing data must be \"0\".
\n") if (@ARGV<2);

my $in = shift;
my $out = shift;
open I, "<", "$in";
open O, ">", "$out";

# Load input file into @gt and %gt.
my $head = <I>;
my (@gt, %gt);
while (<I>) {
	s/\r//;
	chomp;
	my @a = split /\t/;
	push @gt, [@a];
	my $b = shift @a;
	$gt{$b} = \@a;
}
close I;

# Compare every two samples, and load compatible pairs into @line.
my $ns = @gt; # number of sample
my $nl = (@{$gt[0]}-1)/2; # number of loci
my @line; # In the array @line of arrays, every array contains two compatible samples.
foreach my $a (0 .. $ns-2) {
	foreach my $b ($a+1 .. $ns-1) {
		my $d;
		foreach my $c (1 .. $nl) {
			my %na; # number of allele
			$na{$gt[$a][2*$c-1]} = 1;
			$na{$gt[$a][2*$c]} = 1;
			$na{$gt[$b][2*$c-1]} = 1;
			$na{$gt[$b][2*$c]} = 1;
			if (!defined $na{0} and (keys %na) > 2) { # If there are no "0" alleles and number of allele > 2, it is not a compatible locus.
				$d = 1;
				last;
			}
		}
		push @line, [$gt[$a][0], $gt[$b][0]] unless defined $d;
	}
}

# Push certain individual group arrays into @ss, uncertain group arrays into @sn.
my @ss; # samples sure
my @sn; # samples not sure
foreach (@line) {
	if (!defined ${$_}[2]) {
		${$_}[2] = 1; # Mark considered sample pair.
		push my @a, (${$_}[0], ${$_}[1]);
		my @s = @a;
		while (@a) { # The while loop will recursively find $_'s all intercompatible samples, and push them into @s.
			my @c;
			foreach my $v (@a) {
				my @w;
				foreach (@line) { # The foreach loop will find $_'s all compatible samples, and push them into @w.
					if (!defined ${$_}[2]) {
						if (${$_}[0] eq $v or ${$_}[1] eq $v) {
							${$_}[2] = 1;
							push @w, (${$_}[0], ${$_}[1]);
						}
					}
				}
				push @c, @w;		
			}
			@a = @c;
			push @s, @a;
		}
		my %c; # Array @s contain intercompatible samples. The occurrence number of one sample equals the number of its compatible samples. 
		foreach (@s) {
			++$c{$_};
		}
		my %d;
		++$d{$c{$_}} foreach keys %c;
		push @ss, \@s if keys %d == 1; # If the occurrence numbers of all samples are equal, the samples should be one individual.
		push @sn, \@s if keys %d >= 2; # If the numbers are not equal, the samples are uncertain.
	}
}

# Merge samples in @ss into indivuals, and push both indivuals and samples into array @ig.
my @ig; # individual genotype
my $e = 1;
foreach (@ss) {
	my %s;
	foreach (@{$_}) {
		$s{$_} = $gt{$_};
	}
	my @u; # individual genotype
	foreach my $c (0 .. $nl-1) {
		my %t;
		foreach my $b (keys %s) {
			$t{$s{$b}[2*$c]} = 1;
			$t{$s{$b}[2*$c+1]} = 1;
		}
		my @d;# locus genotype
		if (keys %t == 3) { # If there are 3 alleles, push two nonzero alleles into @d.
			foreach (sort {$a <=> $b} keys %t) {
				push @d, $_ if $_ != 0;
			}
		} elsif (keys %t == 2) { # If there are two alleles, ...
			if (defined $t{0}) { # If one of the two is zero, push the nonzero one twice into @d.
				foreach (sort {$a <=> $b} keys %t) {
					push @d, ($_, $_) if $_ != 0;
				}
			} else { # If there is no zero allele, push the two alleles into @d.
				push @d, $_ foreach sort {$a <=> $b} keys %t;
			}
		} elsif (keys %t == 1) { # If there is only one allele, push it twice into @d.
			push @d, (keys %t, keys %t);
		} else {
			warn keys %t;
		}
		push @u, @d;
	}
	unshift @u, sprintf("Indiv%03d", $e);
	++$e;
	my @v;
	push @v, \@u;
	foreach (sort keys %s) {
		my @f;
		push @f, $_;
		push @f, @{$s{$_}};
		push @v, \@f;
	}
	push @ig, \@v;
}

# Print individual groups.
print O "$head\n";
foreach (@ig) {
	foreach (@{$_}) {
		my $a = join "\t", @{$_};
		print O "$a\n";
	}
	print O "\n";
}

# Print uncertain groups.
foreach (@sn) {
	my %s;
	foreach (@{$_}) {
		$s{$_} = $gt{$_};
	}
	print O "Uncertain $e\n";
	++$e;
	foreach (sort keys %s) {
		print O "$_\t";
		print O join "\t", @{$s{$_}};
		print O "\n";
	}
	print O "\n";
}

# Print individuals with only one sample.
foreach (@line) {
	$gt{${$_}[0]} = 0;
	$gt{${$_}[1]} = 0;
}
print O "Individuals With Only One Sample\n";
foreach (sort keys %gt) {
	if ($gt{$_}) {
		print O "$_\t";
		print O join "\t", @{$gt{$_}};
		print O "\n";
	}
}
close O;

