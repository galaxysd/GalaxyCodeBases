#!/usr/bin/perl
use strict;
use warnings;

my $in = shift;
open I, "<", "$in.fa";
open O, ">", "$in.tbl";

my $undefed = "\n\t\t\tnote\tcoding region not determined";

while (<I>) {
	chomp;
	s/>//;
	my @a = split / /;
	my @b = split /-/, $a[0];
	$a[1] =~ s/\[organism=//;
	$a[2] =~ s/\]//;
	my $seq = <I>;
	chomp $seq;
	my $len = length $seq;
	if ($b[0] =~ /PLP/) {
		print O ">Feature $a[0]
<1	>$len	gene
			gene	PLP1
			allele	$a[0]
<1	>61	CDS
			gene	PLP1
			codon_start	3
			product	proteolipid protein 1
62	>$len	intron
			gene	PLP1
";	
	} elsif ($b[0] =~ /16S/) {
		print O ">Feature $a[0]\n";
		print O "<1\t>$len\trRNA\n\t\t\tproduct\tl-rRNA\n\t\t\tnote\t16S ribosomal RNA\n\t\t\tallele\t$b[1]\n";
	} elsif ($b[0] eq 'CYTB') {
		if ($len == 1247) {
			print O ">Feature $a[0]
1	1140	gene
			gene	CYTB
			gene_synonym	cytB
			allele	$b[1]
1	1140	CDS
			gene	CYTB
			gene_synonym	cytB
			codon_start	1
			transl_table	2
			product	cytochrome b
1141	1209	tRNA
			product	tRNA-Thr
1210	>1247	tRNA
			product	tRNA-Pro
";
		} elsif ($len == 1249) {
			print O ">Feature $a[0]
1	1140	gene
			gene	CYTB
			gene_synonym	cytB
			allele	$b[1]
1	1140	CDS
			gene	CYTB
			gene_synonym	cytB
			codon_start	1
			transl_table	2
			product	cytochrome b
1141	1211	tRNA
			product	tRNA-Thr
1212	>1249	tRNA
			product	tRNA-Pro
";
		} elsif ($len == 1233) {
			print O ">Feature $a[0]
<1	1126	gene
			gene	CYTB
			gene_synonym	cytB
			allele	$b[1]
<1	1126	CDS
			gene	CYTB
			gene_synonym	cytB
			codon_start	2
			transl_table	2
			product	cytochrome b
1127	1195	tRNA
			product	tRNA-Thr
1196	>1233	tRNA
			product	tRNA-Pro
";
		} elsif ($len == 853) {
			print O ">Feature $a[0]
1	>853	gene
			gene	CYTB
			gene_synonym	cytB
			allele	$b[1]
1	>853	CDS
			gene	CYTB
			gene_synonym	cytB
			codon_start	1
			transl_table	2
			product	cytochrome b
";
		} else { die '*' x 50,"[x]"; }
	} elsif ( $b[0] eq 'SMCY3' ) {
		print O ">Feature $a[0]
<1	>$len	gene
			gene	SMCY
			allele	$b[1]
<1	>4	CDS
			gene	SMCY
			gene_synonym	KDM5D
			codon_start	2
			product	lysine (K)-specific demethylase 5D
5	>$len	intron
			gene	SMCY
";
	} elsif ( $b[0] eq 'SMCY7_STR_upstream' ) {
		print O ">Feature $a[0]
<1	>$len	gene
			gene	SMCY
			allele	$b[1]
<1	>$len	intron
			gene	SMCY
";
	} elsif ( $b[0] eq 'SMCY7_STR_downstream' ) {
		my $a =  $len - 20;
		my $b = $len - 21;
		print O ">Feature $a[0]
<1	>$len	gene
			gene	SMCY
			allele	$b[1]
<$a	>$len	CDS
			gene	SMCY
			gene_synonym	KDM5D
			codon_start	1
			product	lysine (K)-specific demethylase 5D
<1	$b	intron
			gene	SMCY
";
	} elsif ( $b[0] eq 'DBY7' ) {
		my $a = $len - (278-248);
		my $b = $len - (278-247);
		print O ">Feature $a[0]
<1	>$len	gene
			gene	DBY
			allele	$b[1]
<$a	>$len	CDS
			gene	DBY
			codon_start	1
			product	Y-linked DBY
<1	$b	intron
			gene	DBY
";
	} elsif ( $b[0] eq 'UTY11' ) {
		print O ">Feature $a[0]
<1	>$len	gene
			gene	UTY
			allele	$b[1]
<1	>$len	intron
			gene	UTY
";
	} else {
		print O ">Feature $a[0]\n<1\t>$len\tgene\n\t\t\tgene\t$b[0]\n\t\t\tallele\t$b[1]\n";
		warn "Bad sequence name!";
	}
}
close I;
close O;
my $cmd = "tbl2asn -t Pbe.tpm -i $in.fa -a s2 -V vb";
print "Running:[$cmd]\n";
system($cmd);
#system("ls -l $in.*")
