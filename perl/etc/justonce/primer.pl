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
PRIMER_MIN_SIZE=19
PRIMER_MAX_SIZE=25
PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/usr/local/opt/primer3/share/primer3/primer3_config/
PRIMER_PRODUCT_SIZE_RANGE=30-221
PRIMER_NUM_RETURN=3
=
HEAD

my $in = openfile('u.fa');
while (my $ret = FastaReadNextA($in)) {
	my ($seqname,$genome,$seqdesc) = @$ret;
	#print $seqname," | $seqdesc\n";
	my ($chr,$begin,$end) = split /[:-]/,$seqname;
	my $thePos = $begin +150;
	my $seqlen = length $genome;
	#print "--- $chr,$begin,$end\n";
	print O <<"	ITEM";
SEQUENCE_ID=${seqname}_$thePos
SEQUENCE_TEMPLATE=$genome
SEQUENCE_TARGET=151,1
SEQUENCE_INTERNAL_EXCLUDED_REGION=151,1
SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=101,50,, ; ,,151,50
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
