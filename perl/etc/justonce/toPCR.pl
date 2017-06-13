#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 1.0.0 @ 20120309
=cut
use strict;
use warnings;
use Galaxy::IO;
use Galaxy::IO::FASTA;

my $in = openfile('toPCR.fa.gz');
while (my $ret = FastaReadNext($in)) {
	my ($seqname,$genome) = @$ret;
	print $seqname,"\n";
}

close $in;
