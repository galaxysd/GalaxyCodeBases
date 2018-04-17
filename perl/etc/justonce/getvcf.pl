#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 1.0.0 @ 20180417
=cut
use strict;
use warnings;
use Galaxy::IO;

my $in = openfile('0328.snp.table.16.txt.gz');
chomp(my $t = <$in>);
my @Header = split(/\t/,$t);
print join(", ",@Header),"\n";

while(<$in>) {
	chomp;
	my ($chr,$pos,$Alleles,$Qual,@SampleDat) = split /\t/;
	#print join(", ",$chr,$pos,@SampleDat),"\n";
	my (@SampleGT,@SampleDep);
	for (@SampleDat) {
		my ($TGT,$Deps) = split /;/;
		my @aSampleGT = split /[\/|]/,$TGT;
		my @aSampleDep = split /,/,$Deps;
		my $aSumDep;
		map { $aSumDep += $_ } @aSampleDep;
		push @SampleDep,$aSumDep;
		if (@aSampleGT == 2 and $aSampleGT[0] eq $aSampleGT[1]) {
			push @SampleGT,$aSampleGT[0];
		} else {
			push @SampleGT,join('','x',@aSampleGT);
		}
	}
	print join(", ",$chr,$pos,@SampleGT,@SampleDep),"\n";
}

close $in;
