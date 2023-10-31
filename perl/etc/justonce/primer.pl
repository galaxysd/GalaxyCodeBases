#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 1.0.0 @ 20170612
=cut
use strict;
use warnings;
use Galaxy::IO;
use Galaxy::IO::FASTA qw(FastaReadNextA);

open O,'>','toPCR.p3' or die "Error opening [toPCR.p3]: $!\n";
print O <<'HEAD';
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=24
PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/usr/local/opt/primer3/share/primer3/primer3_config/
PRIMER_PRODUCT_SIZE_RANGE=100-600
PRIMER_MIN_TM=56
PRIMER_OPT_TM=58
PRIMER_MAX_TM=60
PRIMER_NUM_RETURN=1
=
HEAD

my $in = openfile('u.fa');
while (my $ret = FastaReadNextA($in)) {
	my ($seqname,$genome,$seqdesc) = @$ret;
	#print $seqname," | $seqdesc\n";
	#my ($chr,$begin,$end) = split /[:-]/,$seqname;
	#my $thePos = $begin +300;
	my $seqlen = length $genome;
	my $begin = 501;
	my $tglen = $seqlen - 1000;
	#print "--- $chr,$begin,$end\n";
	print O <<"	ITEM";
SEQUENCE_ID=${seqname}_$tglen
SEQUENCE_TEMPLATE=$genome
SEQUENCE_TARGET=501,$tglen
SEQUENCE_INTERNAL_EXCLUDED_REGION=501,$tglen
=
	ITEM
}
print "\n";
close $in;
close O;
print "Run [primer3_core < toPCR.p3 > toPCR.out] now.\n";
__END__
primer3_core -format_output < toPCR.p3.test
primer3_core < toPCR.p3.test

其中一端的引物距离目标位点少于50bp

perl -lane 'print "$F[0]:",$F[1]-300,"-",$F[1]+300' uGRCh37SNPs.txt >uGRCh37SNPs.reg
samtools faidx Homo_sapiens_assembly19.fasta -r uGRCh37SNPs.reg >u.fa &

SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=251,50,, ; ,,301,50
