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
<1\t>$len\tgene
			gene PLP1
<1	62	CDS
			codon_start	3
			product	proteolipid protein 1
			gene PLP1
<1	62	exon
			number	1
			gene PLP1
63	>$len	intron
			number	1
			gene PLP1
			allele	$a[0]
";	
	} elsif ($b[0] =~ /16S/) {
		print O ">Feature $a[0]\n<1\t>$len\tgene\n\t\t\tgene\trrsA\n";
		print O "<1\t>$len\trRNA\n\t\t\tgene\trrsA\n\t\t\tproduct\t16S ribosomal RNA\n\t\t\tallele\t$b[1]\n";
	} elsif ( $b[0] eq 'CytB' ) {
		if ($len == 1247) {
			print O ">Feature $a[0]
1	1140	gene
			gene	CYTB
			allele\t$b[1]
1	1140	CDS
			gene	CYTB
			codon_start	1
			product	cytochrome b
1141	1209	tRNA
			product	tRNA-Thr
1210	>1247	tRNA
			product	tRNA-Pro
";
		} elsif ($len == 853) {
			print O ">Feature $a[0]
1	>853	gene
			gene	CYTB
			allele\t$b[1]
1	>853	CDS
			gene	CYTB
			codon_start	1
			product	cytochrome b
";
		} else { die '*' x 50,"[x]"; }
	} elsif ( $b[0] eq 'SMCY3' ) {
		print O ">Feature $a[0]
<1	>$len	gene
			gene	SMCY
			tallele	$b[1]
54	>$len	intron
			gene	SMCY
			number	4
";
	} elsif ( $b[0] eq 'SMCY7_STR_downstream' ) {
		print O ">Feature $a[0]
<1	>$len	gene
			gene	SMCY7
<1	>$len	intron
			gene	SMCY7
";
	} elsif ( $b[0] eq 'DBY7' ) {
		my $a = $len - (278-248);
		my $b = $len - (278-247);
		print O ">Feature $a[0]
<1	>$len	gene
			gene	DBY
<$a	>$len	CDS
			gene	DBY
			codon_start	1
			number	2
			product	Y-linked DBY
<1	$b	intron
			gene	DBY
";
	} elsif ( $b[0] eq 'DBY7' ) {
		print O ">Feature $a[0]
<1	>$len	gene
			gene	UTY
<99	>$len	intron
			gene	UTY
			number	18
";
	} else {
		print O ">Feature $a[0]\n<1\t>$len\tgene\n\t\t\tgene\t$b[0]\n\t\t\tallele\t$b[1]\n";
	}
}
close I;
close O;
my $cmd = "tbl2asn -t New.tpm -i $in.fa -a s2 -V vb";
print "Running:[$cmd]\n";
system($cmd);
system("ls -l $in.*")

__END__
tbl2asn -t OUT.tpm -i BankIt_submit.fasta -a s2 -V vb -f BankIt_submit.tbl

transl_table	2
note	transl_table=2

cat submit_*.fa > submitall.fa
find submit*.fa|sed 's/\.fa$//'|xargs -n1 perl BankIt_Feature.pl
ls -l *.val



find submit_*.fa|xargs -n1 tbl2asn -t OUT.tpm -a s2 -V vb -i



16S
rRNA            <1..>366
                /product="16S ribosomal RNA"

CytB
gene            1..1140
                /gene="CYTB"
CDS             1..1140
                /gene="CYTB"
                /codon_start=1
                /transl_table=2
                /product="cytochrome b"
                /protein_id="ADK73306.1"
                /db_xref="GI:301340121"
                /translation="MTNIRKSHPLIKIINHSFIDLPAPSNISAWWNFGSLLGVCLILQ
                ILTGLFLAMHYTSDTTTAFSSVTHICRDVNYGWIIRYMHANGASMFFICLYMHVGRGM
                YYGSYTFSETWNIGIMLLFAVMATAFMGYVLPWGQMSFWGATVITNLLSAIPYIGTNL
                VEWIWGGFSVDKATLTRFFAFHFILPFIISALAAVHLLFLHETGSNNPSGITSDSDKI
                PFHPYYTIKDILGLLVLILTLMLLVLFSPDLLGDPDNYIPANPLNTPPHIKPEWYFLF
                AYAILRSIPNKLGGVLALVLSILILAIIPILHTSKQRGMMFRPLSQCLFWLLVADLLT
                LTWIGGQPVEYPFITIGQLASILYFSTLLVLMPISGIIENRLLKW"
tRNA            1141..1209
                /product="tRNA-Thr"
tRNA            complement(1210..>1247)
                /product="tRNA-Pro"


SMCY3
gene            <1..>835
                /gene="SMCY"
intron          54..>835
                /gene="SMCY"
                /number=4



SMCY7_STR_upstream


SMCY7_STR_downstream
gene            <1..>325
                /gene="SMCY7"
intron          <1..>325
                /gene="SMCY7"


DBY7
gene            <1..>278
                /gene="DBY"
CDS             <248..>278
                /gene="DBY"
                /codon_start=2
                /product="Y-linked DBY"
                /protein_id="BAE46786.1"
                /db_xref="GI:77799851"
                /translation="AVKENGRYGRRKQYPISLVLA"
intron          <1..247
                /gene="DBY"



UTY11
gene            <1..>480
                /gene="UTY"
intron          <99..>480
                /gene="UTY"
                /number=18

PLP -> PLP1
gene	<1..>658
	gene	PLP1
CDS             <1..62
                /codon_start=2
                /product="proteolipid protein 1"
                /protein_id="ABN45776.1"
                /db_xref="GI:125522795"
exon            <1..62
                /number=1
intron          63..>658
                /number=1