#!/usr/bin/perl
use strict;
use warnings;

die("
Program: HaploID
Version: 1.0
Release: Sep. 1, 2013\n
Auther: Woody
Consultant: Galaxy\n
Usage: $0 <input> <output>\n
Both input and output files are fasta format.
Use \"N\" for unknown basepair, \"-\" for gap, and \"?\" for missing data.
\n") if (@ARGV<2);

my $in = shift;
my $out = shift;
open I, "<", "$in";

# Read input file.
my %input; # input sequences
my $nchar; # number of characters
my $nin; # number of input sequences
while (<I>) {
	$/ = ">";
	my $seq = <I>;
	chomp $seq;
	$seq =~ s/\s//g;
	$/ = "\n";
	$nchar = length $seq unless $nchar;
	s/>//;
	/^(\S+)\s/;
	push @{$input{$1}}, $seq;
	++$nin;
}
close I;
warn "Read input complete!\n";
warn "No. sequences = $nin\n";

# Make sample consensus sequence.
my %samseq; # sample consensus sequences
my @samid; # sample IDs
foreach my $a (sort keys %input) {
	push @samid, $a;
	my $con = $input{$a}[0];
	my $n = @{$input{$a}} - 1;
	if ($n == 0) {
		$samseq{$a} = $con;
	} else {
		foreach my $b (1 .. $n) {
			for (my $i = 0; $i < $nchar; $i++) {
				my $j = $i + 1;
				if (substr($input{$a}[$b],$i,1) eq "?") {
					next;
				} elsif (substr($input{$a}[$b],$i,1) eq "-") {
					if (substr($con,$i,1) eq "?") {
						substr($con,$i,1) = substr($input{$a}[$b],$i,1);
					} elsif (substr($con,$i,1) eq "-") {
						next;
					} else {
						warn "$a conflicts at $j bp!\n";
					}
				} elsif (substr($input{$a}[$b],$i,1) eq "N") {
					if (substr($con,$i,1) eq "?") {
						substr($con,$i,1) = substr($input{$a}[$b],$i,1);
					} elsif (substr($con,$i,1) eq "-") {
						warn "$a conflicts at $j bp!\n";
					} else {
						next;
					}
				} else {
					if ((substr($con,$i,1) eq "?") or (substr($con,$i,1) eq "N")) {
						substr($con,$i,1) = substr($input{$a}[$b],$i,1);
					} elsif (substr($con,$i,1) eq "-") {
						warn "$a conflicts at $j bp!\n";
					} else {
						if (substr($con,$i,1) eq substr($input{$a}[$b],$i,1)) {
							next;
						} else {
							warn "$a conflicts at $j bp!\n";
						}
					}
				}
			}
		}
		$samseq{$a} = $con;
	}
}
#open O1, ">", "${out}_sample.txt";
#print O1 ">$_\n$samseq{$_}\n" foreach sort keys %samseq;
#close O1;
#warn "Make sample consensus complete and print to ${out}_sample.txt!\n";
my $aa1 = keys %samseq;
warn "No. samples = $aa1\n";

# Compare every two samples, and load compatible pairs into @line.
my @line; # compatible sample pairs
my $nsam = @samid; # total number of samples
foreach my $a (0 .. $nsam-2) {
	foreach my $b ($a+1 .. $nsam-1) {
		my $error;
		for (my $i = 0; $i < $nchar; $i++) {
			if (substr($samseq{$samid[$a]},$i,1) eq "?") {
				next;
			} elsif (substr($samseq{$samid[$a]},$i,1) eq "-") {
				if ((substr($samseq{$samid[$b]},$i,1) eq "?") or (substr($samseq{$samid[$b]},$i,1) eq "-")) {
					next;
				} else {
					$error = 1;
					last;
				}
			} elsif (substr($samseq{$samid[$a]},$i,1) eq "N") {
				if (substr($samseq{$samid[$b]},$i,1) eq "-") {
					$error = 1;
					last;
				} else {
					next;
				}
			} else {
				if ((substr($samseq{$samid[$b]},$i,1) eq "?") or (substr($samseq{$samid[$b]},$i,1) eq "N")) {
					next;
				} elsif (substr($samseq{$samid[$b]},$i,1) eq "-") {
					$error = 1;
					last;
				} else {
					if (substr($samseq{$samid[$a]},$i,1) eq substr($samseq{$samid[$b]},$i,1)) {
						next;
					} else {
						$error = 1;
						last;
					}
				}
			}
		}
		push @line, [$samid[$a], $samid[$b]] unless $error;
	}
}
warn "Compare every two samples complete!\n";
my $aa2 = @line;
warn "No. compatible pairs = $aa2\n";

# Push certain haplotype sample groups into @ss, but one sample many times.
my @unsam; # uncertain sample IDs;
my @ss; # samples sure
my @sn; # samples not sure
$sn[0] = 1;
while (defined $sn[0]) {
	undef @ss;
	undef @sn;
	${$_}[2] = 0 foreach @line;
	foreach (@line) {
		unless (${$_}[2]) {
			${$_}[2] = 1; # Mark considered sample pair.
			push my @a, (${$_}[0], ${$_}[1]);
			my @s = @a;
			while (@a) { # The while loop will recursively find $_'s all intercompatible samples, and push them into @s.
				my @c;
				foreach my $v (@a) {
					my @w;
					foreach (@line) { # The foreach loop will find $_'s all compatible samples, and push them into @w.
						unless (${$_}[2]) {
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
			++$c{$_} foreach @s;
			my %d;
			++$d{$c{$_}} foreach keys %c;
			push @ss, \@s if keys %d == 1; # If the occurrence numbers of all samples are equal, the samples should be one haplotype.
			push @sn, \@s if keys %d >= 2; # If the numbers are not equal, the samples are uncertain.
		}
	}	
	my $bb1 = @ss;
	my $bb2 = @sn;
	warn "No. certain haplotypes = $bb1\n";
	warn "No. uncertain groups  = $bb2\n";
	last unless $bb2;
	foreach (@sn) {
		my %s; # uncertain haplotype sample IDs
		my %sq; # $sq{sample ID} = "?|N" number
		my %nq; # $nl{"?|N" number} = 1
		$s{$_} = 1 foreach @{$_};
		foreach (keys %s) {
			my $a = @{[$samseq{$_} =~ /N|\?/g]};
			$sq{$_} = $a;
			$nq{$a} = 1;
		}
		my $c = ${[sort {$b <=> $a} keys %nq]}[0]; # most "?|N" number in one sample sequence
		foreach my $d (keys %s) {
			if ($sq{$d} == $c) {
				push @unsam, $d;
				foreach my $e (0 .. @line-1) {
					if (defined $line[$e]) {
						if (($line[$e][0] eq $d) or ($line[$e][1] eq $d)) {
							splice @line, $e, 1;
							redo;
						}
					} else {
						last;
					}
				}
			}
		}
		warn "Remove ambiguous samples which have $c \"N|?\".\n";
	}
	my $bb3 = @line;
	warn "No. compatible pairs = $bb3\n";
}
warn "Find haplotype complete!\n";

# Make haplotype consensus.
my @hapseq; # haplotype sequence, @{$hapseq[n][0]} sample IDs, $hapseq[n][1] sequence.
foreach (@ss) {
	my %s; # sample IDs
	foreach (@{$_}) {
		$s{$_} = 1;
	}
	my $con; # consensus sequence;
	my @s; # samples IDs;
	foreach my $a (sort keys %s) {
		$con = $samseq{$a} unless $con;
		for (my $i = 0; $i < $nchar; $i++) {
			if (substr($samseq{$a},$i,1) eq "?") {
				next;
			} elsif ((substr($samseq{$a},$i,1) eq "-") or (substr($samseq{$a},$i,1) eq "N")) {
				if (substr($con,$i,1) eq "?") {
					substr($con,$i,1) = substr($samseq{$a},$i,1);
				} else {
					next;
				}
			} else {
				if ((substr($con,$i,1) eq "?") or (substr($con,$i,1) eq "N")) {
					substr($con,$i,1) = substr($samseq{$a},$i,1);
				} else {
					next;
				}
			}
		}
		push @s, $a;
	}
	push @hapseq, [\@s, $con];
}
warn "Make haplotype consensus complete!\n";

# Print haplotypes with many samples.
open O2, ">", $out;
my $nh = 1; # index of haplotype
foreach (@hapseq) {
	print O2 sprintf(">Haplo%03d\t", $nh), join(" ", @{${$_}[0]}), "\n${$_}[1]\n";
	++$nh;
}

# Print haplotypes which have only one sample and print ambiguous samples.
foreach (@line) {
	delete $samseq{${$_}[0]};
	delete $samseq{${$_}[1]};
}
my $nas = 1; # index of ambiguous sample
my $pas; # print ambiguous samples
foreach (sort @unsam) {
	$pas .= sprintf(">Ambig%03d\t$_\n$samseq{$_}\n", $nas);
	++$nas;
	delete $samseq{$_};
}
foreach (sort keys %samseq) {
	print O2 sprintf(">Haplo%03d\t$_\n$samseq{$_}\n", $nh);
	++$nh;
}
print O2 $pas if $pas;
close O2;
warn "Done and print result to ${out}_haplotype.txt!\n";
