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
PRIMER_PRODUCT_SIZE_RANGE=30-401
PRIMER_NUM_RETURN=1
=
HEAD

my $in = openfile('toPCR.fa.gz');
while (my $ret = FastaReadNextA($in)) {
	my ($seqname,$genome,$seqdesc) = @$ret;
	print $seqname," | $seqdesc\n";
	$seqdesc = (split / /,$seqdesc)[0];
	$seqdesc = (split /:/,$seqdesc)[1];
	my ($begin,$end) = split /-/,$seqdesc;
	my $thePos = (split /_/,$seqname)[1];
	$begin -= $thePos;
	$end -= $thePos;
	print "$thePos, $begin,$end\n";
}

close $in;
close O;
print "Run [primer3_core < toPCR.p3] now.\n";
__END__
primer3_core -format_output < toPCR.p3.test
primer3_core < toPCR.p3.test
