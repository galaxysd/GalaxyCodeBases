#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 1.0.0 @ 20180516
=cut
use strict;
use warnings;

use lib '.';

use Data::Dump qw(ddx);
use Text::NSP::Measures::2D::Fisher::twotailed;
use FGI::GetCPI;
#use Math::BigFloat;

our @Bases;
sub deBayes($) {
	my $p = $_[0];
	my %Dep;
	for my $i (1 .. $#$p) {
		$Dep{$i-1} = $p->[$i];
	}
	my @dKeys = sort { $Dep{$b} <=> $Dep{$a} } keys %Dep;
	if ( $Dep{$dKeys[1]} * 49 > $Dep{$dKeys[0]} ) {	# 2%
		my @rKeys = sort {$a<=>$b} @dKeys[0,1];
		my $gt = join('/',$Bases[$rKeys[0]],$Bases[$rKeys[1]]);
		$p->[0] = $gt;
	}
}

sub getBolsheviks(@) {
	my @dat = map { [split /[;,]/,$_] } @_;
	deBayes($_) for @dat;
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
	#ddx $Bolsheviks,$Hom,\@GTdep,@dat if (keys %GT)>1;
	return [$Bolsheviks,$Hom,\@GTdep];
}

my ($samplesM,$samplesF,$samplesC);

if (@ARGV == 3) {
	($samplesM,$samplesF,$samplesC) = @ARGV;
} elsif (@ARGV == 0) {
	($samplesM,$samplesF,$samplesC) = ('samplesM.tsv','samplesF.tsv','samplesC.tsv');
} else {
	die "Usage: $0 <samplesM> <samplesF> <samplesC>\n";
}

open FM,'<',$samplesM or die "[x]Mom: $samplesM $!\n";
open FF,'<',$samplesF or die "[x]Dad: $samplesF $!\n";
open FC,'<',$samplesC or die "[x]Child: $samplesC $!\n";

my ($logcpi,$lFC,$lFF,$lFM)=(0);
print "# Order: M,F,C\n";

my ($CntA,$CntM) = (0,0);
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
	@Bases = split /,/,$datM[2];
	next if "@tM @tF @tC" =~ /\./;
	#print "@tM\n@datM\n";
	my $retM = getBolsheviks(@tM);
	next unless $retM->[1];
	my $retF = getBolsheviks(@tF);
	#ddx $retM,$retF;
	my @sdatC = map { [split /[;,]/,$_] } @tC;
	my @GTdepC;
	for (@sdatC) {
		for my $i (1 .. $#$_) {
			$GTdepC[$i-1] += $_->[$i];
		}
	}
	my %CntM;
	my @GTM = @{$retM->[2]};
	for my $i (0 .. $#GTM) {
		$CntM{$i} = $GTM[$i];
	}
	my ($x,$y) = sort { $CntM{$b} <=> $CntM{$a} } keys %CntM;
	#ddx [$x,$y,$Bases[$x],$Bases[$y]],\@sdatC,\@GTdepC;
	my $n11 = $retM->[2]->[$x];
	my $n12 = $retM->[2]->[$y];
	my $n21 = $GTdepC[$x];
	my $n22 = $GTdepC[$y];
	next if ($n21+$n22) < 1000;	# skip
	my $GTtC;
	my $twotailedFisher = -1;
	$GTtC = join('/',$Bases[$x],$Bases[$x]);
	my $Cdep = $n21 + $n22;
	#if ($n22 * 199 < $n21) {	# <0.5% = 1:200
	if ($n22/$Cdep < 0.02) {	# minnor < 2%, skip ; depth<10
		next;	# skip
	} else {
		my $n1p = $n11 + $n12;
		my $np1 = $n11 + $n21;
		my $npp = $n1p + $n21 + $n22;
		$twotailedFisher = calculateStatistic(
			n11 => $n11,
			n1p => $n1p,
			np1 => $np1,
			npp => $npp,
		);
		if( (my $errorCode = getErrorCode()) ) {
			die $errorCode, " - ", getErrorMessage();
		} else {
			my ($m,$n) = sort {$a<=>$b} ($x,$y);
			$GTtC = join('/',$Bases[$m],$Bases[$n]);# if $twotailedFisher < 0.05 or $n22 * 49 >= $n21;	# (f0.05 and 0.5%~2%) or >2%
		}
	}
	my $retC = getBolsheviks(@tC);
	my $resM = join(';',$retM->[0],join(',',@{$retM->[2]}));
	my $resF = join(';',$retF->[0],join(',',@{$retF->[2]}));
	my $resC = join(';',$GTtC,join(',',@GTdepC),$twotailedFisher,
						$retC->[0],join(',',@{$retC->[2]})
					);
	my $cret = getcpi(@datM,$resM,$resF,$resC);
	#ddx $cret;
	$logcpi += log($cret->[0]);
	print join("\t",@datM,$resM,$resF,$resC,@$cret,$logcpi/log(10)),"\n";
	++$CntA;
	$CntM = $cret->[1];
}

close FM; close FF; close FC;

my $lgcpi = sprintf("%f",$logcpi/log(10));
my @cpiPart = split /\./,$lgcpi;
my ($main,$texp) = (1,$cpiPart[0]);
if (scalar @cpiPart == 2) {
	my $mainstr;
	if ($logcpi < 0) {
		$mainstr = "-0.$cpiPart[1]";
	} else {
		$mainstr = "0.$cpiPart[1]";
	}
	$main = 10*exp($mainstr*log(10));
	$texp -= 1;
}

print "# Total: $CntA, Mis: $CntM, CPI: ",join('',$main,'e',$texp),"\n";

__END__

Order M,F,C

Canceled:
+放弃贝叶斯结果，2.5%以上就是杂合。双亲只保留多数结果。
x子代单个样品深度<1000的，整行扔掉。
x子代，both >0.5% and chi^2<0.05，才算杂合。

All the words below is provided by the client, original text:

1. 先读母亲和儿子的文件，只选map质量大于100，并且母亲是醇和snp的位点
我:
我:
2， 如果是多次测序，少数服从多数，并吧多数的加起来。把三个人的genotype定住
3. 确定基因型：（计算母亲与儿子snp的差异的卡方的pvalue，要求孩子的alt rate大于0.5%，并且卡方pvalue<0.05,则取与母亲不同的genotype，否则，取和母亲一致的genotype）

4，计算CPI
注意事项;1,母亲的genotype只用纯和的（以vcf结果为准），2，孩子的深度必须大于1000X
结果： snp位点  孩子genotype(合起来的#ref/合起来的#alt) 母亲genotype 父亲genotype
