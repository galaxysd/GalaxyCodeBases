#!/bin/env perl
use strict;
use warnings;
#use Data::Dump qw(dump ddx);

unless (@ARGV){
	print "perl $0 <soaps.lst> <outprefix(.stat & .ins)>\n";	# soaps.nfo can store the file size of soap. Too late, useless.
	exit;
}

my ($fqlst,$statout) = @ARGV;
my (%DATrbrf,%nfo,%insD);

open LST,'<',$fqlst or die "[x]Error opening $fqlst: $!\n";
while (<LST>) {
	chomp;
	my ($PESE,$sample,$lib,$FL,$ReadLen,$nfofpath)=split /\t/;
	$nfo{$sample}{$lib}{$FL}=[$nfofpath,$PESE,$ReadLen];
	$DATrbrf{Lane}{$sample}{$lib}{$FL}={ALL => [0,0,0,0,0,{},{},{},{},0], Summary => [0,0,0], };
	# "ReadsOut, BPOut, MisSum, TrimedReads, TrimedBP, misMatchReads, Reads@Hit, BP@Hit, IndelReads, BadLines"
	# "Total_Reads(=Total_Pairs*2), Paired, Singled/Alignment"
	$DATrbrf{Lib}{$sample}{$lib}={ALL => [0,0,0,0,0,{},{},{},{},0], Summary => [0,0,0], } unless exists $DATrbrf{Lib}{$sample}{$lib};	# may be faster ?
	$DATrbrf{Sample}{$sample}={ALL => [0,0,0,0,0,{},{},{},{},0], Summary => [0,0,0], } unless exists $DATrbrf{Sample}{$sample};
	$insD{$sample}{$lib}{$FL}=['',''];
}

sub combineC($) {
	my $href=$_[0];
	if ($href and %$href) {
		my (@str,$m);
		$m = (sort {$a<=>$b} keys %$href)[-1];
		for (1..$m) {
			push @str,$$href{$_}||0;
		}
		return \join('|',@str);
	} else {return \'.';}
}
sub combineLineA($) {
	my $ref=$_[0];
	return join ',',@{$$ref{Summary}},@{$$ref{ALL}}[0..4],${&combineC($$ref{ALL}[5])},${&combineC($$ref{ALL}[6])},${&combineC($$ref{ALL}[7])},${&combineC($$ref{ALL}[8])},$$ref{ALL}[9];
}
sub combineLine($) {
	my $ref=$_[0];
	return join ',',@$ref[0..4],${&combineC($$ref[5])},${&combineC($$ref[6])},${&combineC($$ref[7])},${&combineC($$ref[8])};
}
sub sumcsv ($$) {
	my ($href,$stref)=@_;
	return if $$stref eq '.';
	my %new = map {split /:/} split(/,/,$$stref);
	$$href{$_} += $new{$_} for (keys %new);
}
sub sumup ($$$) {
	my ($PESE,$sum,$item)=@_;
	my ($a,$m)=($$sum{Summary},$$sum{ALL});
	my ($b,$n)=($$item{Summary},$$item{ALL});
	if ($PESE eq 'PE') {
		$$a[0] += $$b[0]+$$b[0];
		$$a[1] += $$b[1];
		$$a[2] += $$b[2];
	} else {
		$$a[0] += $$b[0];
		$$a[2] += $$b[1];
	}
	$$m[-1] += $$n[-1];	# BadLines
	for my $chr (keys %$item) {
		next if $chr eq 'Summary';
		#$$sum{$_}=[0,0,0,0,{},{},{},{}] unless $$sum{$_};
		($m,$n)=($$sum{$chr},$$item{$chr});
		$$sum{$chr}=$m=[0,0,0,0,0,{},{},{},{}] unless $m;
		$$m[$_] += $$n[$_] for (0..4);
		&sumcsv($$m[$_],\$$n[$_]) for (5..8);
#warn "$chr\n";
	}
}

my %NFO;
for my $sample (sort keys %nfo) {
	for my $lib (keys %{$nfo{$sample}}) {
		for my $FL (keys %{$nfo{$sample}{$lib}}) {
			my ($nfofpath,$PESE,$ReadLen)=@{$nfo{$sample}{$lib}{$FL}};
			%NFO=();
			#print "[$_]\n";
			open NFO,'<',"$nfofpath" or (warn "[!]Error opening $nfofpath: $!\n" and next);
warn "[!]$nfofpath\n";
			while (<NFO>) {
				next if /^(#|$)/;
				chomp;
				my ($key,@values)=split /\t/;
				$key='ALL' if $key eq '_ALL_';	# for old format
				$NFO{$key}=\@values;
			}
			close NFO;
			&sumup($PESE,$DATrbrf{Sample}{$sample},\%NFO);
			&sumup($PESE,$DATrbrf{Lib}{$sample}{$lib},\%NFO);
			&sumup($PESE,$DATrbrf{Lane}{$sample}{$lib}{$FL},\%NFO);
			$insD{$sample}{$lib}{$FL}=[@{$NFO{Summary}}[3,4]] if $PESE eq 'PE';
		}
	}
}
#ddx \%DATrbrf;
my %Chr;
open O,'>',$statout.'.stat' or die "[x]Error opening $statout.stat: $!\n";
print O "#Summary\tSubItemOrder: Total_Reads,Aligned_Pairs,Aligned_Single,ReadsOut,BPOut,MisSum,TrimedReads,TrimedBP,misMatchReads|ASC,Reads\@Hit|ASC,BP\@Hit|ASC,IndelReads|ASC,BadLines\n";
open INS,'>',$statout.'.ins' or die "[x]Error opening $statout.ins: $!\n";
print INS "#Sample\tLib\tLane\tMode(p%),Lsd,Rsd,InsAvg,STD\tInsMin,InsMax\n";
my ($Rsample,$Rlib,$RFL);
for my $sample (sort keys %nfo) {
	$Rsample=&combineLineA($DATrbrf{Sample}{$sample});#join ',',@{$DATrbrf{Sample}{$sample}{Summary}},@{$DATrbrf{Sample}{$sample}{ALL}}[0..3],${&combineC($DATrbrf{Sample}{$sample}{ALL}[4])},${&combineC($DATrbrf{Sample}{$sample}{ALL}[5])},${&combineC($DATrbrf{Sample}{$sample}{ALL}[6])},${&combineC($DATrbrf{Sample}{$sample}{ALL}[7])},$DATrbrf{Sample}{$sample}{ALL}[8];
	for (keys %{$DATrbrf{Sample}{$sample}}) {
		next if /(Summary|ALL)/;
		++$Chr{$_};
	}
	for my $lib (sort keys %{$nfo{$sample}}) {
		$Rlib=&combineLineA($DATrbrf{Lib}{$sample}{$lib});
		for my $FL (sort keys %{$nfo{$sample}{$lib}}) {
			$RFL=&combineLineA($DATrbrf{Lane}{$sample}{$lib}{$FL});
			print O join("\t",'ALL',$sample,$Rsample,$lib,$Rlib,$FL,$RFL),"\n";
			print INS join("\t",$sample,$lib,$FL,@{$insD{$sample}{$lib}{$FL}}),"\n";
		}
	}
}
close INS;
print O "\n#ByChr\tSubItemOrder: ReadsOut,BPOut,TrimedReads,TrimedBP,misMatchReads|ASC,Reads\@Hit|ASC,BP\@Hit|ASC,IndelReads|ASC\n";
for my $Chr (sort keys %Chr) {
	for my $sample (sort keys %nfo) {
		$Rsample=&combineLine($DATrbrf{Sample}{$sample}{$Chr});
		for my $lib (sort keys %{$nfo{$sample}}) {
			$Rlib=&combineLine($DATrbrf{Lib}{$sample}{$lib}{$Chr});
			for my $FL (sort keys %{$nfo{$sample}{$lib}}) {
				$RFL=&combineLine($DATrbrf{Lane}{$sample}{$lib}{$FL}{$Chr});
				my $ChrID=$Chr;
				$ChrID =~ s/^;//;
				print O join("\t",$ChrID,$sample,$Rsample,$lib,$Rlib,$FL,$RFL),"\n";
			}
		}
	}
}
close O;
warn "[!] Done !\n";

__END__
对./to9311/2soap/soaps.stat：
上半部分：
#Summary        SubItemOrder: Total_Reads,Aligned_Pairs,Aligned_Single,ReadsOut,BPOut,TrimedReads,TrimedBP,misMatchReads|ASC,Reads@Hit|ASC,BP@Hit|ASC,IndelReads|ASC,BadLines
按个体统计：grep ^ALL to9311/2soap/soaps.stat | awk '{print $2"\t"$3}' |sort|uniq
按文库统计：grep ^ALL to9311/2soap/soaps.stat | awk '{print $4"\t"$5}' |sort|uniq
按lane统计：grep ^ALL to9311/2soap/soaps.stat | awk '{print $6"\t"$7}'

下半部分：
#ByChr  SubItemOrder: ReadsOut,BPOut,TrimedReads,TrimedBP,misMatchReads|ASC,Reads@Hit|ASC,BP@Hit|ASC,IndelReads|ASC
按个体统计：grep ^Chr01 to9311/2soap/soaps.stat | awk '{print $2"\t"$3}' |sort|uniq
按文库统计：grep ^Chr01 to9311/2soap/soaps.stat | awk '{print $4"\t"$5}' |sort|uniq
按lane统计：grep ^Chr01 to9311/2soap/soaps.stat | awk '{print $6"\t"$7}'

以Chr01为例，其他染色体替换。或用“Chr”提取全部。

misMatchReads|ASC示例：
“300347|184857|0|222”，mismatch=1的有300347，mismatch=2的有184857，没有mismatch=3的，mismatch=4的有222。
IndelReads|ASC 同上。

BP@Hit|ASC与Reads@Hit|ASC类似，但只分为1..9和>=10共计10种类别。
根据抽取一个soap结果统计结果，不同hit数的分布是类似双曲线趋势下降的。
