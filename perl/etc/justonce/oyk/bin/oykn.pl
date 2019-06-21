#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 1.0.0 @ 20180516
=cut
use strict;
use warnings;
use POSIX;
use FindBin qw($RealBin);
if ($FindBin::VERSION < 1.51) {
	warn "[!]Your Perl is too old, thus there can only be ONE `bsuit` file in your PATH. [FindBin Version: $FindBin::VERSION < 1.51]\n\n"
}
FindBin::again();
use lib "$RealBin/../";
require FGI::rsCPI;

#use Data::Dump qw(ddx);
my $DBsuffix = '../db/nippt7274.tsv';

my @Modes = qw(CHIP PCR);
my %Mode = map { $_ => 1 } @Modes;
my @Parentages = qw(DUO TRIO);
my %Parentage = map { $_ => 1 } @Parentages;
my $Verbose = 0;

die "Usage: $0 <Chip|PCR> <Duo|Trio> <mother> <father> <child> <outprefix>

[Mode]: PCR mode require 2 repeats while Chip mode require 1.
[Parentage]: A duo test involves the child and an alleged father. 例如代孕。
In contrast, a trio test involves the mother, child, and the alleged father.\n" if @ARGV<6;
$Verbose = 1 if @ARGV == 7;
my $theMode = uc shift;
unless (exists $Mode{$theMode}) {
	die "[x]Mode can only be:[",join(',',@Modes),"].\n";
}
my $theParentage = uc shift;
unless (exists $Parentage{$theParentage}) {
	die "[x]Parentage can only be:[",join(',',@Parentages),"].\n";
}

sub getTPE(@) {
	my @AFs = @_;
	my ($p1,$p2)=(0,0);
	for my $i (0 .. $#AFs) {
		my $p1i = $AFs[$i]*(1-$AFs[$i])*(1-$AFs[$i]);
		$p1 += $p1i;
		for my $j (($i+1) .. $#AFs) {
			my $p2i = $AFs[$i]*$AFs[$i]*$AFs[$j]*$AFs[$j]*(4-3*$AFs[$i]-3*$AFs[$j]);
			$p2 += $p2i;
		}
	}
	return $p1-$p2;
}
sub getDPE(@) {
	my @AFs = @_;
	my ($p1,$p2)=(0,0);
	for my $i (0 .. $#AFs) {
		my $p1i = $AFs[$i]*$AFs[$i]*(1-$AFs[$i])*(1-$AFs[$i]);
		$p1 += $p1i;
		for my $j (($i+1) .. $#AFs) {
			my $tb = 1-$AFs[$i]-$AFs[$j];
			my $p2i = 2*$AFs[$i]*$AFs[$j]*$tb*$tb;
			$p2 += $p2i;
		}
	}
	return $p1+$p2;
}
our (%Markers,%MarkerAF);
open DB,'<',"$RealBin/$DBsuffix" or die $?;
while (<DB>) {
	next if /^#/;
	chomp;
	my ($rs,$chr,$pos,$chr38,$pos38,@d) = split /\t/,$_;
	$MarkerAF{$rs} = {@d};
	my @AFs = values %{$MarkerAF{$rs}};
	#my $sum = eval join '+',@AFs;
	#print "$rs: $sum\n" if $sum != 1;
	my $TPE = getTPE(@AFs);
	my $DPE = getDPE(@AFs);
	$Markers{$rs} = [$chr,$pos,$DPE,$TPE];
}
#ddx \%Markers;
#ddx \%MarkerAF;
close DB;
#print getDPE(0.5,0.5),"\n";print getTPE(0.5,0.5),"\n";
#0.125
#0.1875

our @Bases;
sub deBayes($) {
	my $p = $_[0];
	my %Dep;
	for my $i (1 .. $#$p) {
		$Dep{$i-1} = $p->[$i];
	}
	#ddx %Dep;
	my @dKeys = sort { $Dep{$b} <=> $Dep{$a} } keys %Dep;
	if ( @dKeys>1 and $Dep{$dKeys[1]} >= $Dep{$dKeys[0]} * 0.02) {	# 2%
		my @rKeys = sort {$a<=>$b} @dKeys[0,1];
		my $gt = join('/',$Bases[$rKeys[0]],$Bases[$rKeys[1]]);
		$p->[0] = $gt;
	} elsif (@dKeys>1 && ($Dep{$dKeys[1]} < $Dep{$dKeys[0]} * 0.02) && ($Dep{$dKeys[1]} > $Dep{$dKeys[0]} * 0.001)){
		$p->[0] = "NA";
	} elsif (@dKeys == 1 or ($Dep{$dKeys[1]} <= $Dep{$dKeys[0]} * 0.001)){
		my $gt = join('/',$Bases[$dKeys[0]],$Bases[$dKeys[0]]);
		$p->[0] = $gt;
	}
}
sub deBayes2($) {
	my $p = $_[0];
	my %Dep;
	for my $i (1 .. $#$p) {
		$Dep{$i-1} = $p->[$i];
	}
	#ddx %Dep;
	my @dKeys = sort { $Dep{$b} <=> $Dep{$a} } keys %Dep;
	if ( @dKeys>1 and $Dep{$dKeys[1]}  >= $Dep{$dKeys[0]} * 0.1 ) {	# 10%
		my @rKeys = sort {$a<=>$b} @dKeys[0,1];
		my $gt = join('/',$Bases[$rKeys[0]],$Bases[$rKeys[1]]);
		$p->[0] = $gt;
	} elsif (@dKeys>1 && ($Dep{$dKeys[1]}  > $Dep{$dKeys[0]} * 0.01) && ($Dep{$dKeys[1]}  < $Dep{$dKeys[0]} * 0.1)){
		$p->[0] = "NA";
	}
}

sub getBolsheviks(@) {
	my $type = shift;
	my @dat = map { [split /[;,]/,$_] } @_;
	if ($type) {
		deBayes($_) for @dat;
	} else {
		deBayes2($_) for @dat;
	}
	#ddx \@dat;

	my (%GT);
	for (@dat) {
		++$GT{$_->[0]};
	}
	if (defined $GT{NA}){
		return ["NA",0,"NA"];
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

sub getequal(@) {
	my $type = shift;
	my @dat = map { [split /[;,]/,$_] } @_;
#ddx \@dat;
	if ($type) {
		deBayes($_) for @dat;
	} else {
		deBayes2($_) for @dat;
	}
	my (%GT);
	for (@dat) {
		++$GT{$_->[0]};
	}
	for (values %GT) {
		return 1 if $_ == 2;	# 两个样品GT一致
	}
	return 0;
}
sub getrio(@) {
	my @dat = map { [split /[;,]/,$_] } @_;
	#ddx \@dat;
	my @ret;
	for (@dat) {
		shift @$_;
		my @depth = @$_;
		my @sortAsc = sort {$a <=> $b} @depth;
		#ddx \@depth;
		my ($sum,$cnt)=(0,0);
		for (@depth) {
			$sum += $_;
			++$cnt if $_;
		}
		if ($cnt == 3) {
			push @ret,[1,$sortAsc[0]/$sum,$sortAsc[1]/$sum,join(',',@depth)];
			$ret[-1]->[-1] = $ret[-1]->[-1] . ',T' if $Verbose;
		} elsif ($cnt == 4) {
			push @ret,[1,($sortAsc[0]+$sortAsc[1])/$sum,$sortAsc[2]/$sum,join(',',@depth)];
			$ret[-1]->[-1] = $ret[-1]->[-1] . ',Q' if $Verbose;
		} elsif ($cnt == 1) {
			push @ret,[0,0,0,join(',',@depth)];
			$ret[-1]->[-1] = $ret[-1]->[-1] . ',S' if $Verbose;
		} elsif ($cnt == 2) {
			push @ret,[0,0,0,join(',',@depth)];
			$ret[-1]->[-1] = $ret[-1]->[-1] . ',D' if $Verbose;
		} else {
			die;
		}
	}
	my ($n,$x,$xx,$y,$yy,$depstr,@depstrs)=(0,0,0,0,0,'');
	if (@ret == 0) {
		die;
	} else {
		for (@ret) {
			my $flag = shift @$_;
			if ($flag) {
				$x += $$_[0];
				$xx += $$_[0]*$$_[0];
				$y += $$_[1];
				$yy += $$_[1]*$$_[1];
				++$n;
			}
			push @depstrs,$$_[2]; 
		}
		$depstr = join(';',@depstrs);
	}
	#ddx \@ret;
	return ($n,$x,$xx,$y,$yy,$depstr);
}
sub tstat(%) {
	my %d = @_;
	unless ($d{'n'}) {
		return ('NA','NA');
	}
	my $mean1 = $d{'x'}/$d{'n'};
	my $mean2 = $d{'y'}/$d{'n'};
	my $std1 = sqrt($d{'xx'}/$d{'n'} - $mean1*$mean1);
	my $std2 = sqrt($d{'yy'}/$d{'n'} - $mean2*$mean2);
	my $srt1 = join(' ± ',$mean1,$std1);
	my $srt2 = join(' ± ',$mean2,$std2);
	return ($srt1,$srt2);
}
sub printExp($) {
	my $lnV = $_[0]/log(10);
	my $lnInt = floor($lnV);
	my $lnExt = $lnV - $lnInt;
	my $prefix = exp($lnExt*log(10));
	my $str = join('e',$prefix,$lnInt);
	return $str;
}

my $mother=shift;
my $father=shift;
my $child=shift;
my $outprefix=shift;
open FM,'<',$mother or die "[x]Mom: $!\n";
open FF,'<',$father or die "[x]Dad: $!\n";
open FC,'<',$child or die "[x]Child: $!\n";

open OC,'>',"$outprefix.cpie" or die "[x]$outprefix.cpie: $!\n";
open OT,'>',"$outprefix.trio" or die "[x]$outprefix.trio: $!\n";
open OR,'>',"$outprefix.tHM" or die "[x]$outprefix.tsv: $!\n";

my ($logcpi,$spe,$trioN,$lFC,$lFF,$lFM)=(0,0,0);
my (%trioM,%trioF,%trioC);
for (qw(n x xx y yy)) {
	$trioM{$_} = 0;
	$trioF{$_} = 0;
	$trioC{$_} = 0;
}
print "# Order: M,F,C\n";

while (<FM>) {
	chomp;
	chomp($lFC = <FC>);
	chomp($lFF = <FF>);
	my @datM = split /\t/;
	my @datF = split /\t/,$lFF;
	my @datC = split /\t/,$lFC;
	#my ($chr,undef,$bases,$qual,@data) = split /\t/;
	next if $datM[3] !~ /\d/ or $datM[3] < 100;
	next if $datF[3] !~ /\d/ or $datF[3] < 100;
	next if $datC[3] !~ /\d/ or $datC[3] < 100;
	die if $datM[0] ne $datC[0] or $datF[0] ne $datC[0];
	my @tM = splice @datM,4;
	my @tF = splice @datF,4;
	my @tC = splice @datC,4;
	@Bases = split /,/,$datM[2];	# $bases = ref,alt
	next if $Bases[1] eq '.';
	next if "@tM @tF @tC" =~ /\./;

	my $retM = getBolsheviks(0,@tM);
	my $retF = getBolsheviks(0,@tF);
	#ddx $retM if $retM->[1];
	my (@rM,@rF,@rC);
	if (@Bases > 2) {
		#ddx \@datM,\@datF,\@datC;
# oyka.pl:220: (
#   ["SNP5501", 501, "C,A,G", 3763.92],
#   ["SNP5501", 501, "C,A,G", 3763.92],
#   ["SNP5501", 501, "C,A,G", 3763.92],
# )
		#ddx \@tM,\@tF,\@tC;
# oyka.pl:221: (
#   ["C/A;28,26,0", "C/A;28,26,0"],
#   ["C/G;35,0,37", "C/G;35,0,37"],
#   ["C/A;36,48,0", "C/A;36,48,0"],
# )
		@rM = getrio(@tM);
		@rF = getrio(@tF);
		@rC = getrio(@tC);
		#ddx \@rM,\@rF,\@rC;
#   [
#     2,
#     0.560606060606061,
#     0.157139577594123,
#     0.606060606060606,
#     0.183654729109275,
#     "37,55,40,T;37,55,40,T",
#   ],
		if ($rM[0]+$rF[0]+$rC[0] >0) {
			for (qw(n x xx y yy)) {
				$trioM{$_} += shift @rM;
				$trioF{$_} += shift @rF;
				$trioC{$_} += shift @rC;
			}
#   ["37,55,40,T;37,55,40,T"]
			#ddx (\%trioF,\%trioM,\%trioC);
			++$trioN;
			if ($retM->[1]) {
				print OR join("\t",@datM[0,2],$rM[0],$rF[0],$rC[0]),"\n";
				my $str = join('=',$retM->[0],join(',',@{$retM->[2]}));
				$rM[0] = join('.',$rM[0],$str,'HM');
			}
			if ($retF->[1]) {
				my $str = join('=',$retF->[0],join(',',@{$retF->[2]}));
				$rF[0] = join('.',$rF[0],$str,'HF');
			}
			print OT join("\t",$trioN,@datM[0,2],$rM[0],$rF[0],$rC[0]),"\n";
		}
	}
	my $check_dep = 1;
	for (@tM){
		my @info = split /[;,]/,$_;
		my $sum;
		for my $i(1..scalar @info - 1){
			$sum += $info[$i];
		}
		if ($sum > 50){
			$check_dep *= 1;
		}else{
			$check_dep *= 0;
		}
	}
	for (@tF){
		my @info = split /[;,]/,$_;
		my $sum;
		for my $i(1..scalar @info - 1){
			$sum += $info[$i];
		}
		if ($sum > 50){
			$check_dep *= 1;
		}else{
			$check_dep *= 0;
		}
	}
	next if ($check_dep == 0);

	#T/T;6,2245      C/C;1698,0
	#print "> @tM , @tF , @tC\n@datM\n";
	#my $retM = getBolsheviks(0,@tM);
	next unless $retM->[1];
	#my $retF = getBolsheviks(0,@tF);
	next if ($retM->[0] eq "NA" or $retF->[0] eq "NA");
	#ddx $retM,$retF;
	my $xx = getequal(0,@tM);
	my $yy = getequal(0,@tF);
	my $zz = getequal(1,@tC);
	my $t=$xx*$yy*$zz;
	my $REP = 1;
	if ($theMode eq 'PCR') {
		next unless $t;
	} else {
		$REP = 2 if $t;
	}
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
	my ($n12,$n22);
	my $n11 = $retM->[2]->[$x];
	my $n21 = $GTdepC[$x];
	if (defined $y) {
		$n12 = $retM->[2]->[$y];
		$n22 = $GTdepC[$y];
	} else {
		$n12 = 0;
		$n22 = 0;
	}
	next unless defined $n22;
	if ($theMode eq 'PCR') {
		next if ($n21+$n22) < 200;
	} elsif ($theMode eq 'CHIP') {
		next if ($n21+$n22) < (100 * $REP);
	}
	my $GTtC;
	$GTtC = join('/',$Bases[$x],$Bases[$x]);
	my $Cdep = $n21 + $n22;

	my $retC = getBolsheviks(1,@tC);
	#ddx $retM,$retF,$retC;
	next if ($retC->[0] eq "NA");
	my @fgeno=split /\//,$retF->[0];
	my @mgeno=split /\//,$retM->[0];
	my @cgeno=split /\//,$retC->[0];

	if ($theMode eq 'PCR') {
		next if $fgeno[0] eq $fgeno[1] and $mgeno[0] eq $mgeno[1] and $mgeno[0] eq $fgeno[0] and (($retM->[2][0]>0 and $retM->[2][1]>0) or ($retF->[2][0]>0 and $retF->[2][1]>0));
	}

	my @mnum=@{$retM->[2]};
	my @fnum=@{$retF->[2]};	
	my $resM = join(';',$retM->[0],join(',',@{$retM->[2]}));
	my $resF = join(';',$retF->[0],join(',',@{$retF->[2]}));
	my $resC = join(';',$retC->[0],join(',',@{$retC->[2]}));
#	my $resC = join(';',$GTtC,join(',',@GTdepC),$twotailedFisher,
#						$retC->[0],join(',',@{$retC->[2]})
#					);
	my $cret;
	if ($theParentage eq 'TRIO') {
		if ($theMode eq 'CHIP') {
			next if $mgeno[0] eq $mgeno[1] && $cgeno[0] eq $cgeno[1] && $mgeno[0] eq $cgeno[0];
		}
		$cret = getcpiT(@datM,$resM,$resF,$resC);
	} elsif ($theParentage eq 'DUO') {
		my $mGT = $mgeno[0];
		if ($cgeno[0] eq $cgeno[1]) {
			my $d1 = $retC->[2]->[0];
			my $d2 = $retC->[2]->[1];
			my $flag = 0;
			$flag = 1 if $d1 > 500;
			$flag = 1 if $d2 > 500;
			next if $flag == 0;
		}
		$cret = getcpiD(@datM,$resM,$resF,$resC);
	}
# @datM,$resM,$resF,$resC: (
#   "SNP4593",
#   501,
#   "G,C",
#   13170,
#   "C/C;0,176",
#   "C/C;0,206",
#   "G/C;26,314",
# )
	#ddx $cret;
	$logcpi += log($cret->[0]);
	$spe += log(1-$cret->[1]);
	print OC join("\t",@datM,$resM,$resF,$resC,@$cret,$logcpi/log(10),$spe/log(10)),"\n";
}

close FM; close FF; close FC;

my $sCPI = printExp($logcpi);
my $sPCPE = printExp($spe);

print OC "# CPI: $sCPI\n";
print OC "# CPE: 1-$sPCPE\n";

print "CPI: $sCPI\n";
print "CPE: 1-$sPCPE\n";

if ($trioN) {
	my @stM = tstat(%trioM);
	my @stF = tstat(%trioF);
	my @stC = tstat(%trioC);
	#ddx (\@stM,\@stF,\@stC);
	print OT "# M1: $stM[0] , M2: $stM[1]\n";
	print OT "# F1: $stF[0] , F2: $stF[1]\n";
	print OT "# C1: $stC[0] , C2: $stC[1]\n";

	print "M1: $stM[0] , M2: $stM[1]\n";
	print "F1: $stF[0] , F2: $stF[1]\n";
	print "C1: $stC[0] , C2: $stC[1]\n";
}

close OC;close OT;close OR;

__END__
grep '[ACTG],[ATCG],[ATCG]' *.tsv|grep '[1-9],[1-9],[1-9]'
./oykn.pl chip Trio s385M1.tsv s385F1.tsv s385C.tsv ss

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
