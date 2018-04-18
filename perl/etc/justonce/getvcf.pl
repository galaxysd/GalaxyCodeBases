#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 1.0.0 @ 20180417
=cut
use strict;
use warnings;
use Galaxy::IO;
use Data::Dump qw(ddx);

my $in = openfile('0328.snp.table.16.txt.gz');
chomp(my $t = <$in>);
my @Header = split(/\t/,$t);
my @SampleIDs = splice @Header,4;
my $qGroupCnt = @SampleIDs/4;
die "[x]Sample count $qGroupCnt is not 4n !\n" if $qGroupCnt != int($qGroupCnt);
print join(", ","$qGroupCnt: ",@SampleIDs),"\n";

my (%Result,%GroupTGT);

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
	for my $i (1 .. $qGroupCnt) {
		my $flag = 1;
		my %gTGTs;
		for my $j (4*($i-1) .. (4*$i-1)) {
			if ($SampleDep[$j] < 15) {	# min depth is 15
				$flag = 0;
				last;
			} else {
				++$gTGTs{$SampleGT[$j]}
			}
		}
		my $cp = "$chr\t$pos";
		if ($flag == 0) {
			$GroupTGT{$cp}{$i} = '.';	# less depth
			next;
		}
		my @TGTs = keys(%gTGTs);
		if (@TGTs == 1 and length($TGTs[0])!=1) {
			$GroupTGT{$cp}{$i} = '-';	# 4 being same, but Het
			next;
		} else {
			++$Result{$i}{@TGTs};
			++$Result{$i}{-1};
			if (@TGTs == 1) {
				$GroupTGT{$cp}{$i} = $TGTs[0];
			}
		}
	}
	#print join(", ",$chr,$pos,@SampleGT,@SampleDep),"\n";
	ddx \%Result;
}
close $in;

ddx \%GroupTGT;
ddx \%Result;
