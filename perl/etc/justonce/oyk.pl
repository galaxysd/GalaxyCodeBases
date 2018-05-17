#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 1.0.0 @ 20180516
=cut
use strict;
use warnings;
use Galaxy::IO;
use Data::Dump qw(ddx);

use Text::NSP::Measures::2D::Fisher::twotailed;

my $n11 = 1968;
my $n12 = 4;
my $n21 = 1975;
my $n22 = 15;

my $n1p = $n11 + $n12;
my $np1 = $n11 + $n21;
my $npp = $n1p + $n21 + $n22;

my $twotailed_value = calculateStatistic(
	n11 => $n11,
	n1p => $n1p,
	np1 => $np1,
	npp => $npp,
);

if( (my $errorCode = getErrorCode()) ) {
	print STDERR $errorCode, " - ", getErrorMessage();
} else {
	print getStatisticName, "value for bigram is ", $twotailed_value, "\n";
}

sub getBolsheviks(@) {
	my @dat = map { [split /[;,]/,$_] } @_;
	my (%GT);
	for (@dat) {
		++$GT{$_->[0]};
	}
	my $Bolsheviks = (sort {$GT{$b} <=> $GT{$a}} keys %GT)[0];
	my @GTdep;
	for (@dat) {
		next if $_->[0] ne $Bolsheviks;
		for my $i (1 .. $#$_) {
			$GTdep[$i-1] += $_->[$i];
		}
	}
	my @GTs = split /[\/|]/,$Bolsheviks;
	my $Hom = 0;
	$Hom = 1 if $GTs[0] eq $GTs[1];
	if ((keys %GT)>1) {
		ddx $Bolsheviks,$Hom,\@GTdep,@dat;
	}
	return [$Bolsheviks,$Hom,\@GTdep];
}

open FM,'<','samplesM.tsv' or die "[x]Mom: $!\n";
open FF,'<','samplesF.tsv' or die "[x]Dad: $!\n";
open FC,'<','samplesC.tsv' or die "[x]Child: $!\n";

my ($lFC,$lFF,$lFM);

while (<FM>) {
	chomp;
	chomp($lFC = <FC>);
	chomp($lFF = <FF>);
	my @datM = split /\t/;
	my @datF = split /\t/,$lFF;
	my @datC = split /\t/,$lFC;
	#my ($chr,undef,$bases,$qual,@data) = split /\t/;
	next if $datM[3] < 100;
	next if $datF[3] < 100;
	next if $datC[3] < 100;
	die if $datM[0] ne $datC[0];
	my @tM = splice @datM,4;
	my @tF = splice @datF,4;
	my @tC = splice @datC,4;
	next if "@tM @tF @tC" =~ /\./;
	#print "@tM\n@datM\n";
	my $retM = getBolsheviks(@tM);
	next unless $retM->[1];
	my $retF = getBolsheviks(@tF);
	ddx $retM,$retF;
}

close FM; close FF; close FC;

__END__

Order M,F,C

All the words below is provided by the client, original text:

1. 先读母亲和儿子的文件，只选map质量大于100，并且母亲是醇和snp的位点
我:
我:
2， 如果是多次测序，少数服从多数，并吧多数的加起来。把三个人的genotype定住
3. 确定基因型：（计算母亲与儿子snp的差异的卡方的pvalue，要求孩子的alt rate大于0.5%，并且卡方pvalue<0.05,则取与母亲不同的genotype，否则，取和母亲一致的genotype）

4，计算CPI
注意事项;1,母亲的genotype只用纯和的（以vcf结果为准），2，孩子的深度必须大于1000X
结果： snp位点  孩子genotype(合起来的#ref/合起来的#alt) 母亲genotype 父亲genotype
