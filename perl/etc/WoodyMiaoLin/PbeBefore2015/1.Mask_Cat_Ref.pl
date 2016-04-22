#!/usr/bin/perl
use strict;
use warnings;

open REF, "<", "/bak/seqdata/genomes/Felis_catus-6.2/Felis_catus-6.2.fa";
open GENE, "-|", "gzip -dc /share/users/miaolin/5.Leopard_cat/f.Gronau_2011_analysis/3.Data_Filters/4.UCSC_Track_Filters/ensGene.txt.gz";
open EXON, "-|", "gzip -dc /share/users/miaolin/5.Leopard_cat/f.Gronau_2011_analysis/3.Data_Filters/4.UCSC_Track_Filters/ensemblSource.txt.gz";
open OUT, "|-", "gzip -9c >Felis_catus-6.2_masked.fa.gz";

my %chr = (
	chrA1 => ["gi|362110686|gb|CM001378.1|", 239302903],
	chrA2 => ["gi|362110683|gb|CM001379.1|", 169043629],
	chrA3 => ["gi|362110681|gb|CM001380.1|", 142459683],
	chrB1 => ["gi|362110674|gb|CM001381.1|", 205241052],
	chrB2 => ["gi|362110670|gb|CM001382.1|", 154261789],
	chrB3 => ["gi|362110654|gb|CM001383.1|", 148491654],
	chrB4 => ["gi|362110649|gb|CM001384.1|", 144259557],
	chrC1 => ["gi|362110644|gb|CM001385.1|", 221441202],
	chrC2 => ["gi|362110642|gb|CM001386.1|", 157659299],
	chrD1 => ["gi|362110640|gb|CM001387.1|", 116869131],
	chrD2 => ["gi|362110638|gb|CM001388.1|", 89822065],
	chrD3 => ["gi|362110636|gb|CM001389.1|", 95741729],
	chrD4 => ["gi|362110634|gb|CM001390.1|", 96020406],
	chrE1 => ["gi|362110632|gb|CM001391.1|", 63002102],
	chrE2 => ["gi|362110629|gb|CM001392.1|", 64039838],
	chrE3 => ["gi|362110625|gb|CM001393.1|", 43024555],
	chrF1 => ["gi|362110622|gb|CM001394.1|", 68669167],
	chrF2 => ["gi|362110619|gb|CM001395.1|", 82763536],
	chrX => ["gi|362110616|gb|CM001396.1|", 126427096],
);

my %refseq;
while (<REF>) {
	chomp;
	/Felis catus breed Abyssinian chromosome (\w+)/;
	my $name = "chr$1";
	$/ = ">";
	my $seq = <REF>;
	chomp $seq;
	$seq =~ s/\s//g;
	my $len = length $seq;
	warn "$len\t$chr{$name}[1]" if $len != $chr{$name}[1];
	$/ = "\n";
	$refseq{$name} = $seq;
}
close REF;
warn "Read sequences complete!\n";

my %exon;
while (<EXON>) {
	chomp;
	my @a = split /\t/;
	$exon{$a[0]} = 1 if $a[1] eq "protein_coding";
}
close EXON;
warn "Read EXON complete!\n";

while (<GENE>) {
	my @a = split /\t/;
	next unless grep(/^$a[2]$/, keys %chr);
	next unless $exon{$a[1]};
	my @sta = split /,/, $a[9]; # Starts of exons
	my @end = split /,/, $a[10]; # Ends of exons
	warn $_ if @sta != @end;
	my $i = @sta - 2; # Number of exons minus one
	foreach (0 .. $i) {
		warn $_ if $sta[$_] >= $end[$_];
		my $s = $sta[$_]; # Start of filter
		my $e = $end[$_]; # End of filter
		foreach ($s .. $e) {
			substr($refseq{$a[2]}, $_, 1) = lc substr($refseq{$a[2]}, $_, 1);
		}
	}
}
close GENE;
warn "Mask GENE complete!\n";

warn "#Chr\tNo. Total\tNo. Masked\tNo. Remained\n";
foreach my $u (sort keys %refseq) {
	print OUT ">$u $chr{$u}[0]\n";
	my $i;
	my $i0;
	my $i1;
	foreach my $v (0 .. $chr{$u}[1]-1) {
		++$i;
		if (substr($refseq{$u}, $v, 1) =~ /[A-Z]/) {
			++$i0;
			print OUT "N";
		} elsif (substr($refseq{$u}, $v, 1) =~ /[a-z]/) {
			++$i1;
			print OUT uc(substr($refseq{$u}, $v, 1));
		} else {
			warn;
		}
		if ($i == 100) {
			print OUT "\n";
			$i = 0;
		}
	}
	print OUT "\n";
	warn "$u\t$chr{$u}[1]\t$i0\t$i1\n";
}
close OUT;
warn "Done!\n";
