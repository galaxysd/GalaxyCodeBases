#!/usr/bin/perl
use strict;
use warnings;

open I, "-|", "samtools view /share/users/huxs/work/catbtail/bam/FCAP072_H2LGFCCXX_L3.aln.bam |awk '{print \$9}'";

my $num;
my $sum;
while (<I>) {
	chomp;
	if ($_ > 0 and $_ < 1000) {
		$sum += $_;
		$num += 1;
	}
}

close I;
print $sum, "\t", $num, "\t", $sum/$num, "\n";
