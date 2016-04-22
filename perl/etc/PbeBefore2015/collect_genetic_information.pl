#!/usr/bin/perl
use strict;
use warnings;

open I0, "<", "../a.information/Haplotype_Region.txt";
my %hap_reg;
while (<I0>) {
	chomp;
	my @c = split /\t/;
	next if !$c[0];
	$hap_reg{$c[0]} = $c[1];
}
close I0;

open I1, "<", "../a.information/Sample_MT_X_Y.txt";
my %sample_gt; # Sample, MT, X, Y, ATP8, Autosome.
my %sample_re; # Sample, MT, X, Y, ATP8, Autosome.
while (<I1>) {
	chomp;
	my @c = split /\t/;
	my $d = shift @c;
	$sample_gt{$d} = [@c];
	$sample_re{$d} = [$hap_reg{$c[0]}, $hap_reg{$c[1]}, $hap_reg{$c[2]}];
}
close I1;

open I2, "<", "/share/users/miaolin/5.Leopard_cat/3.PBE_paper_2015/g.ATP8_haplotype_tree/b.modified_alignment/f.modified_alignment/ATP8_haplotype_location_sampleID.txt";
while (<I2>) {
	chomp;
	next if !/^Pbe/;
	my @c = split /\t/;
	my @d = split / /, $c[2];
	foreach (@d) {
		push @{$sample_gt{$_}}, $c[0];
		push @{$sample_re{$_}}, $hap_reg{$c[0]};
	}
}
close I2;
foreach (keys %sample_gt) {
	if (@{$sample_gt{$_}} == 3) {
		push @{$sample_gt{$_}}, 0;
	} elsif (@{$sample_gt{$_}} == 4) {
		next;
	} else {
		warn @{$sample_gt{$_}};
	}
}

foreach (keys %sample_re) {
	if (@{$sample_re{$_}} == 3) {
		push @{$sample_re{$_}}, 0;
	} elsif (@{$sample_re{$_}} == 4) {
		next;
	} else {
		warn @{$sample_re{$_}};
	}
}

open I3, "<", "../a.information/Sample_AutoRegion.txt";
while (<I3>) {
	chomp;
	my @c = split /\t/;
	push @{$sample_gt{$c[0]}}, $c[1];
	push @{$sample_re{$c[0]}}, $c[2];
}
close I3;
foreach (keys %sample_gt) {
	if (@{$sample_gt{$_}} == 4) {
		push @{$sample_gt{$_}}, 0;
	} elsif (@{$sample_gt{$_}} == 5) {
		next;
	} else {
		warn @{$sample_gt{$_}};
	}
}

foreach (keys %sample_re) {
	if (@{$sample_re{$_}} == 4) {
		push @{$sample_re{$_}}, 0;
	} elsif (@{$sample_re{$_}} == 5) {
		next;
	} else {
		warn @{$sample_re{$_}};
	}
}

foreach (sort keys %sample_gt) {
	foreach (@{$sample_gt{$_}}) {
		$_ = "n/a" unless $_;
	}
	my %region;
	foreach (@{$sample_re{$_}}) {
		$_ = "n/a" unless $_;
		$region{$_} = 1;
	}
	delete $region{"n/a"};
	if ((keys %region) == 1) {
		push @{$sample_re{$_}}, (keys %region);
	} elsif ((keys %region) == 2) {
		if ($region{Northern}) {
			if ($region{Iriomote}) {
				push @{$sample_re{$_}}, "Iriomote";
			} elsif ($region{"China/Indochina"}) {
				push @{$sample_re{$_}}, "China/Indochina";
			} elsif ($region{Taiwan}) {
				push @{$sample_re{$_}}, "Taiwan";
			} elsif ($region{Amur}) {
				push @{$sample_re{$_}}, "Amur";
			} elsif ($region{Sunda}) {
				push @{$sample_re{$_}}, "Hybrid";
			} elsif ($region{Philippines}) {
				push @{$sample_re{$_}}, "Hybrid";
			} elsif ($region{Southern}) {
				push @{$sample_re{$_}}, "Hybrid";
			} else {
				push @{$sample_re{$_}}, "CHECK!";
			}
		} elsif ($region{Southern}) {
			if ($region{Iriomote}) {
				push @{$sample_re{$_}}, "Hybrid";
			} elsif ($region{"China/Indochina"}) {
				push @{$sample_re{$_}}, "Hybrid";
			} elsif ($region{Taiwan}) {
				push @{$sample_re{$_}}, "Hybrid";
			} elsif ($region{Amur}) {
				push @{$sample_re{$_}}, "Hybrid";
			} elsif ($region{Sunda}) {
				push @{$sample_re{$_}}, "Sunda";
			} elsif ($region{Philippines}) {
				push @{$sample_re{$_}}, "Philippines";
			} else {
				push @{$sample_re{$_}}, "CHECK!";
			}
		}
	} else {
		if ($region{Northern} and $region{Southern}) {
			push @{$sample_re{$_}}, "Hybrid";
		} else {
			push @{$sample_re{$_}}, "CHECK!";
		}
	}
}

open O, ">", "sample_genotype_region.txt";
print O "#SampleID\tgtMT\tgtX\tgtY\tgtATP8\tqAutosome(K=2)\tregionMT\tregionX\tregionY\tregionATP8\tregionAutosome\tGeneticRegion\n";
foreach (sort keys %sample_gt) {
	print O $_, "\t";
	print O join("\t", @{$sample_gt{$_}}), "\t";
	print O join("\t", @{$sample_re{$_}}), "\n";
}
close O;
