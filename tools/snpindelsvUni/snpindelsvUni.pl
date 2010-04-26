#!/bin/env perl
use strict;
use warnings;
use lib '/share/raid010/resequencing/soft/lib';
use Time::HiRes qw ( gettimeofday tv_interval );
use File::Basename;
use Galaxy::ShowHelp;

$main::VERSION=0.0.1;

###
my $U_m=1;	# CDS区影响权重pre bp
my $U_s=1.1;	# 平均移码影响值（pre CDS）
#my $U_b=0.12;	# 标准化调控区比例（含启动子、增强子加权值）
my $U_n=0.9;	# 调控区影响权重pre bp

my $SNP_HET_Ratio=1.6;		# 杂合SNP系数
my $Indel_HET_Ratio=0.9;	# 杂合Indel系数
my $SV_Mixture_Ratio=0.75;	# 杂合Indel系数
my %SV_Weights=(			# 复合sv拆分系数
	Transposion => 0.5,
);							# 转座的权重算1/2
my $SV_Confidence=0.25;		# SV可靠度 / sv不确定度:v

my $SNP_r=1000;
my $Indel_r=1000;
my $SV_r=10;
###

our $opts='i:r:g:p:n:s:o:c:vb';
our($opt_i, $opt_r, $opt_g, $opt_p, $opt_n, $opt_s, $opt_v, $opt_b, $opt_o, $opt_c);

our $desc='';
our $help=<<EOH;
\t-i Genome Info [ChrID\\tLen]
\t-r RegEx for ChrID (chromosome_\\d+)
\t-g GFF file for CDS length
\t-c CFG file for SNP, Indel & SV
\t-p SNP file
\t-n Indel file
\t-s SV file
\t-o Output file ({CFG}.out)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_r = 'chromosome_\d+' unless $opt_r;
$opt_o = $opt_c.'.out' unless $opt_o;

if ($opt_c) {
	if (-s $opt_c) {
		warn "[!]Reading from [$opt_c].\n";
		open C,'<',$opt_c or warn "[x]Error: $!\n";
		while (<C>) {
			chomp;
			my ($k,$v)=split /\t/;
			$opt_p=$v if uc($k) eq 'SNP';
			$opt_n=$v if uc($k) eq 'INDEL';
			$opt_s=$v if uc($k) eq 'SV';
		}
		close C;
	}
}
print STDERR "From [$opt_i][$opt_g][$opt_p][$opt_n][$opt_s] to [$opt_o] with to [$opt_r]\n";
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
$|=1;
open( IN,'<',$opt_i) or die "[x]Error: $!\n";
my $GenomeLen;
while (<IN>) {
	chomp;
	my ($chrid,$len)=split /\t/;
	next unless /$opt_r/;
	$GenomeLen += $len;
	print "$chrid,$len\t$GenomeLen\n" if $opt_v;
}
close IN;
print "[!]GenomeLen: $GenomeLen\n";

open( IN,'<',$opt_g) or die "[x]Error: $!\n";
my ($CDSLen,$GeneLen,$CDSCount,$GeneCount);
while (<IN>) {
	chomp;
	my ($seqname, $source, $primary, $start, $end, $score, $strand, $frame, $groups) = split /\t/;
	next unless $seqname =~ /$opt_r/;
	$primary = uc $primary;
	if ($primary eq 'CDS') {
		$CDSLen += $end-$start+1;
		++$CDSCount;
	} elsif ($primary eq 'GENE') {
		$GeneLen += $end-$start+1;
		++$GeneCount;
	}
}
close IN;
my $U_a=$CDSLen/$GeneLen;	# CDS比例;
my $U_b=(750*$GeneCount+4*$CDSCount)/$GeneLen;	# 标准化调控区比例（含启动子、增强子加权值）
print "[!]GeneLen: $GeneLen, GeneCount: $GeneCount, CDSLen: $CDSLen, CDSCount: $CDSCount
[!]a=$U_a, b=$U_b\n";

=pod
chromosome_10	1226	A	M	36	C	28	4	5	A	33	3	6	11	0.0285714	1.54545	0	166
1.	Chromosome name
2.	Position of locus
3.	Nucleotide at corresponding locus of reference sequence
4.	Genotype of sequencing sample
5.	Quality value
6.	nucleotide with the highest probability(first nucleotide)
7.	Quality value of the nucleotide with the highest probability
8.	Number of supported reads that can only be aligned to this locus
9.	Number of all supported reads that can be aligned to this locus
10.	Nucleotide with higher probability
11.	Quality value of nucleotide with higher probability
12.	Number of supported reads that can only be aligned to this locus
13.	Number of all supported reads that can be aligned to this locus
14.	Total number of reads that can be aligned to this locus
15.	RST p-value
16.	Estimated copy number for this locus
17.	Presence of this locus in the dbSNP database. 1 refers to presence and 0 refers to inexistence
18.	The distance between this locus and another closest SNP
=cut
open( IN,'<',$opt_p) or die "[x]Error: $!\n";
my ($CountSNP,$i,$U_SNP);
while (<IN>) {
	chomp;
	my ($chr, $pos, $ref, $snp, $Q, $best, undef, undef, undef, $second, undef, undef, undef,
		$depth, $rst, $CN) = split /\t/;
	next unless $chr =~ /$opt_r/;
	if ($snp =~ /[ATCG]/i) {
		$i=1;
	} else {
		$i=$SNP_HET_Ratio;
	}
	$CountSNP += $i/$CN;
}
close IN;
$U_SNP=$SNP_r*$CountSNP*($U_a*$U_m+$U_b*$U_n)/$GeneLen;
print "[!]SNP Value: $CountSNP -> $U_SNP\n";

=pod
chromosome_1	2936	I2	TC	*	homo	27	15	17
1.	Chromosome
2.	Position
3.	Indel type and number
4.	Bases
5.	Strand (+, positive; -, negative; *, both positive and negative)
6.	Homozygosis or heterozygosis
7.	Average quality
8.	The number of supported read pair
9.	The number of all crossed read pair
=cut
open( IN,'<',$opt_n) or die "[x]Error: $!\n";
my ($CountIndel,$len,$U_Indel);
while (<IN>) {
	chomp;
	my ($chr, $pos, $IDv, $bases, $strand, $HomHet, $Q) = split /\t/;
	next unless $chr =~ /$opt_r/;
	$len=length $bases;
	$len *= $Indel_HET_Ratio if $HomHet eq 'hete';
	$CountIndel += $len;
}
close IN;
$U_Indel=$Indel_r*$CountIndel*(2*$U_s + $U_a*$U_m + 3*$U_b*$U_n)/(3*$GeneLen);
print "[!]Indel Value: $CountIndel -> $U_Indel\n";

=pod
chromosome_10	Deletion	7455	8413	7461	8386	925	0-0-0	FR	4
1.	Chromosome name
2.	Type of structure variation
3.	Minimal value of start position in cluster
4.	Maximal value of end position in cluster
5.	Estimated start position of this structure variation
6.	Estimated end position of this structure variation
7.	Length of SV
8.	Breakpoint of SV (only for insertion)
9.	Unusual matching mode (F refers to align with forward sequence, R refers to align with reverse sequence)
10.	number of paired-end read which support this structure variation
=cut
open( IN,'<',$opt_s) or die "[x]Error: $!\n";
my (%CountSV,@types,$CountSVA,$U_SV);
while (<IN>) {
	chomp;
	my ($chr, $type, $left, $right, $start, $end, $len) = split /\t/;
	next unless $chr =~ /$opt_r/;
	unless ($type =~ /\+/) {
		$CountSV{$type} += $len;
	} else {
		@types = split /\+/,$type;
		$CountSV{$_} += $len*$SV_Mixture_Ratio for @types;
	}
}
close IN;
print "[!]SV Details: ";
for my $type (sort keys %CountSV) {
	print "$type: $CountSV{$type}   ";
	$CountSV{$type} *=$SV_Weights{$type} if exists $SV_Weights{$type};
	$CountSVA += $CountSV{$type} * $SV_Confidence;
}
$U_SV=$SV_r*$CountSVA*(2*$U_s + $U_a*$U_m + 3*$U_b*$U_n)/(3*$GeneLen);
print "\n[!]SV Value: $CountSVA -> $U_SV\n";

my $U_mark=$U_SNP+$U_Indel+$U_SV;
print '[!]Summary: ',$U_mark," = $U_SNP+$U_Indel+$U_SV\n";
open O,'>',$opt_o or die "[x]Error: $!\n";
print O "SNP:\t$U_SNP\nIndel:\t$U_Indel\nSV:\t$U_SV\nSummary:\t$U_mark\n";
close O;

#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
#print "\nPlease use the following command to batch qsub:\033[32;1m
#find ${opt_o}/ -name '*.sh' | while read ll; do qsub -l vf=2G -cwd \$ll; done\n\033[0;0m\n";
__END__
./snpindelsvUni.pl -i sbi1.genome -g Sbi1.4.gff3 -p Ji_snp.txt -n Ji.indel -s Ji_sv.txt -v
