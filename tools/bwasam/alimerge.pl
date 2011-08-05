#!/bin/env perl
use strict;
use warnings;
use Data::Dump qw(dump ddx);

unless (@ARGV){
	print "perl $0 <outprefix(.stat & .ins)> <stator info files>\n";	# soaps.nfo can store the file size of soap. Too late, useless.
	exit;
}

my $statout = shift @ARGV;
my (%DATrbrf,%nfo,%insD);

my $Dathref={ALL => [0,0,0,0,0,{},{},{},{},0], Summary => [0,0,0], };
# "ReadsOut, BPOut, MisSum, TrimedReads, TrimedBP, misMatchReads, Reads@Hit, BP@Hit, IndelReads, BadLines"
# "Total_Reads(=Total_Pairs*2), readsPaired, Singled/Alignment"
my @InsStr;

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
sub combineJ($) {
	my $href=$_[0];
	if ($href and %$href) {
		my @str;
		for (sort {$a<=>$b} keys %$href) {
			push @str,join(':',$_,$$href{$_});
		}
		return \join(',',@str);
	} else {return \'.';}
}
sub combineLineA($) {
	my $ref=$_[0];
	return join "\t",@{$$ref{Summary}},@{$$ref{ALL}}[0..4],${&combineC($$ref{ALL}[5])},${&combineC($$ref{ALL}[6])},${&combineC($$ref{ALL}[7])},${&combineJ($$ref{ALL}[8])},$$ref{ALL}[9];
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
my $PESE='PE';
open INS,'>',$statout.'.ins' or die "[x]Error opening $statout.ins: $!\n";
while(my $nfofpath=shift @ARGV) {
    next unless -f $nfofpath;
	%NFO=();
	print STDERR "Read [$nfofpath]\n";
	print INS "#[$nfofpath]\n";
	open NFO,'<',"$nfofpath" or (warn "[!]Error opening $nfofpath: $!\n" and next);
	while (<NFO>) {
		next if /^(#|$)/;
		chomp;
		my ($key,@values)=split /\t/;
		$key='ALL' if $key eq '_ALL_';	# for old format
		$NFO{$key}=\@values;
	}
	close NFO;
	&sumup($PESE,$Dathref,\%NFO);
	push @InsStr,${$NFO{Summary}}[3] if $PESE eq 'PE';
	#ddx \%NFO;
}
#ddx \%DATrbrf;
ddx $Dathref;
ddx \@InsStr;

print INS "\n#Mode(p%),Lsd,Rsd,InsAvg,STD\n";
print INS "$_\n" for @InsStr;
close INS;
#__END__
my %Chr;
open O,'>',$statout.'.stat' or die "[x]Error opening $statout.stat: $!\n";
print O "#Summary:\n#Total_Reads	Aligned_Pairs,Aligned_Single\tReadsOut\tBPOut\tMisSum\tTrimedReads\tTrimedBP\tmisMatchReads|ASC\tReads\@Hit|ASC\tBP\@Hit|ASC\tIndelReads|ASC\tBadLines\n";
my $RFL=&combineLineA($Dathref);
print O $RFL,"\n";
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
