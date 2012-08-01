#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Data::Dump qw(ddx);
use Galaxy::IO;
use Galaxy::SeqTools;

die "Usage: $0 <bcgv bcf> <out>\n" if @ARGV<1;
my $bcfs=shift;
my $outfs=shift;

open OP,'>',$outfs.'.tped' or die "Error opening $outfs.tped : $!\n";
open OM,'>',$outfs.'.MinorAllele' or die "Error opening $outfs.MinorAllele : $!\n";
$t = "# In: [$bcfs], Out: [$outfs]\n";
print O $t;
print $t;





print "Prepare $outfs.tfam [and $outfs.phe]. And then:\np-link --tfile $outfs --reference-allele $outfs.MinorAllele --fisher --out $outfs.p --pheno $outfs.phe --all-pheno --model --cell 0\n";
__END__
grep -hv \# radseq.gt > radseq.tfam
