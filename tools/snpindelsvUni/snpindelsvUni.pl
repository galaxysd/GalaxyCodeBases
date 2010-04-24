#!/bin/env perl
use strict;
use warnings;
use lib '/share/raid010/resequencing/soft/lib';
use Time::HiRes qw ( gettimeofday tv_interval );
use File::Basename;
use Galaxy::ShowHelp;

$main::VERSION=0.0.1;

###
my $SNP_HET_Ratio=1.6;		# 杂合SNP系数
my $Indel_HET_Ratio=0.9;	# 杂合Indel系数
###

our $opts='i:r:g:p:n:s:o:vb';
our($opt_i, $opt_r, $opt_g, $opt_p, $opt_n, $opt_s, $opt_v, $opt_b, $opt_o);

our $desc='';
our $help=<<EOH;
\t-i Genome Info [ChrID\\tLen]
\t-r RegEx for ChrID (chromosome_\\d+)
\t-g GFF file for CDS length
\t-p SNP file
\t-n Indel file
\t-s SV file
\t-o Output file
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_r = 'chromosome_\d+' unless $opt_r;

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
my ($CDSLen,$GeneLen);
while (<IN>) {
	chomp;
	my ($seqname, $source, $primary, $start, $end, $score, $strand, $frame, $groups) = split /\t/;
	next unless $seqname =~ /$opt_r/;
	$primary = uc $primary;
	$CDSLen += $end-$start+1 if $primary eq 'CDS';
	$GeneLen += $end-$start+1 if $primary eq 'GENE';
}
close IN;
print "[!]GeneLen: $GeneLen, CDSLen: $CDSLen\n";

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
my ($CountSNP,$i);
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
print "[!]SNP Value: $CountSNP\n";

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
my ($CountIndel,$len);
while (<IN>) {
	chomp;
	my ($chr, $pos, $IDv, $bases, $strand, $HomHet, $Q) = split /\t/;
	next unless $chr =~ /$opt_r/;
	$len=length $bases;
	$len *= $Indel_HET_Ratio if $HomHet eq 'hete';
	$CountIndel += $len;
}
close IN;
print "[!]Indel Value: $CountIndel\n";

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
my ($CountSV,$len);
while (<IN>) {
	chomp;
	my ($chr, $type, $left, $right, $start, $end, $len) = split /\t/;
	next unless $chr =~ /$opt_r/;
	$len *= $Indel_HET_Ratio if $HomHet eq 'hete';
	$CountIndel += $len;
}
close IN;
print "[!]Indel Value: $CountIndel\n";


#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
#print "\nPlease use the following command to batch qsub:\033[32;1m
#find ${opt_o}/ -name '*.sh' | while read ll; do qsub -l vf=2G -cwd \$ll; done\n\033[0;0m\n";
__END__
./snpindelsvUni.pl -i sbi1.genome -g Sbi1.4.gff3 -p Ji_snp.txt -n Ji.indel -s Ji_sv.txt -v
