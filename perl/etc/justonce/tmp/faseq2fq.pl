#!/opt/gentoo/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Data::Dump qw (ddx);

my $infq = $ARGV[0];
my $infa = $ARGV[1];
my $outf = $ARGV[2];

print "FQ.gz:[$infq] <- FA:[$infa] to Out.fq.gz:[$outf]\n";

#$infq = '/share/users/liuyc/Project_Tiger_Population_Genomics/3_Tiger_Data_Align/PSMC/pti096_diploid.fq.gz';
#$infa = '/share/users/liuyc/Project_Tiger_Population_Genomics/3_Tiger_Data_Align/Tiger_SNP_fa/pti096_GQ20.fa';
#$infa = 't.fa';

my %Genome;

open FA,'<',$infa;
$/ = '>';
while(<FA>) {
	chomp;
	next if $_ eq '';
	#print "[$_]\n";
	my @dat = split /\n/;
	my $id = shift @dat;
	my $seq = join('',@dat);
	#print "[$id][$seq]\n";
	#die;
	$Genome{$id}=$seq;
}
close FA;

$/ = "\n";

my $seqin = Bio::SeqIO->new(
	-file   => "gzip -dc $infq |",
	-format => 'fastq',
	);

my $seqout = Bio::SeqIO->new(
	-file   => "| gzip -9c >$outf",
	-format => 'fastq',
	);	

while (my $inseq = $seqin->next_seq) {
	my $id = $inseq->{'display_id'};
	if (not exists $Genome{$id}) {
		warn "Err $id\n";
		ddx $inseq;
	}
	$inseq->{'seq'} = $Genome{$id};
	$seqout->write_seq($inseq);
}

__END__
./faseq2fq.pl /share/users/liuyc/Project_Tiger_Population_Genomics/3_Tiger_Data_Align/PSMC/diploid_fq/pti096_diploid.fq.gz t.fa out.fq.gz &

ddx $inseq;
->{}
	alphabet       => "dna",
	display_id     => "scaffold1001",
	namespace      => "sanger",
	seq            => "actgagtcacccaggtgcccctgaattacacattctaaaggggttAAAATGGTGAATTTTGTGTCAGACACATTTAGCATACTAAAAGTCACACATTCTTCTGTATCTACAAGAAAAG",
	_mapping       => [1, 1],
	_meta          => {
						DEFAULT => [
										30,
	                                    32,
...skipping...
										60,
										33,
									],
						trace   => [],
						},
	_nowarnonempty => undef,
	_ranges        => undef,
	_root_verbose  => 0,
	_seq_length    => undef,
