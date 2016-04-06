#!/usr/bin/env perl
use strict;
use warnings;

die "Usage: $0 <SimLst> <AnalyseRes> [Shift=10]\n" if @ARGV <2;

my $total=shift;
my $result=shift;
my $shift=shift;
$shift = 10 unless defined $shift;
open TT,$total or die $!;
open RT,$result or die $!;
my %hash;

while(<TT>){
	chop;
	my @a=split;
	#print $a[3]."\n";
	next unless(/\w/);
	for( my $kk=$a[3]-$shift;$kk<$a[3]+$shift;$kk++){
		$hash{$kk}=1;
	}

}
close TT;

sub docheck($$$) {
	my ($chr,$start,$end) = @_;
	my $ret = 0;
	($start,$end) = sort {$a <=> $b} ($start,$end);
	for my $p ($start .. $end) {
		if (exists $hash{$p}) {
			$ret = 1;
			last;
		}
	}
	return $ret;
}

while(<RT>){
	chomp;
	my @a=split;
	my ($cnt,$chr,$start,$end,$chr2,$start2,$end2)=(1);
	if (@a == 11) { # invers & tandem_dup
		($chr,$start,$end) = @a[5,6,7];
	} elsif (@a == 14 or @a == 15) { # ins & del_ins
		($chr,$start,$end) = @a[5,6,7];
		($chr2,$start2,$end2) = @a[9,10,11];
		$cnt=2 if $chr2 ne '-';
	} elsif (@a == 13) { # transl_inter
		($chr,$start,$end) = @a[5,6,6];
		($chr2,$start2,$end2) = @a[8,9,9];
		$cnt=2 if $chr2 ne '-';
	} else {
		die "Unknown format:[$_]\n"
	}
	#print $a[2]."\n";
	my $ret = docheck($chr,$start,$end);
	$ret = docheck($chr2,$start2,$end2) if ($cnt==2) and ($ret==0);
	print "$_\n" if $ret;
}
close RT;

__END__
grep \> simout_*.Ref.fa | sed 's/^simout_m//'|sed 's/.Ref.fa:>/\t/'|sed 's/Ref_/Ref:/g'|sed 's/Vir_/Vir:/'|sed 's/R_/R:/'|sed 's/_/ /g'|cat -n >simed.lst

$ ./overlapMeerkat.pl simed.lst all.variants|awk '{print $1"\t"$2}'|uniq|sort|uniq
11	invers
11	invers_f
11	invers_r
11	tandem_dup
13	transl_inter
14	del_ins
15	del_inso
15	inso
15	inss

wc -l all.variants	117	# cat Meerkat/*2*.variants > all.variants
./overlapMeerkat.pl simed.lst all.variants 100 |wc -l	82

$ wc -l simed.lst
     700 simed.lst

$ wc -l Meerkat/*2*.variants
       8 Meerkat/bwa2gx.variants
      86 Meerkat/bwa2hg18.variants
      23 Meerkat/meth2gx.variants
     117 total

./overlapMeerkat.pl simed.lst Meerkat/bwa2gx.variants |wc -l	3
./overlapMeerkat.pl simed.lst Meerkat/meth2gx.variants |wc -l	5
./overlapMeerkat.pl simed.lst Meerkat/bwa2hg18.variants |wc -l	71

__DAT__
invers_f	NA	273_0	31	1	chr18	11929626	11930843	1217	1000	BP:TE_TE
invers_r	NA	270	32	1	chr18	11929606	11930823	1217	1000	BP:TE_TE
tandem_dup	NA	50	97	2	chr18	12616569	12617117	548	-47	BP:_TE

del_ins	FoSTeS	1574_0	5	1	chr18	15924715	18626457	2701741	-	-	-	41	BP:_TE
invers_f	NA	1189	10	1	chr18	61614378	61615078	700	1000	BP:_
invers_r	NA	1358	7	1	chr18	61614358	61615058	700	1000	BP:_
tandem_dup	NA	772	47	1	chr18	22015346	22016298	952	20	BP:TE_TE
tandem_dup	NA	300	101	2	chr18	43138397	43138712	315	-49	BP:_

del_inso	FoSTeS	815_0/474_0	5/18	20/2	chr6	147769523	147769543	19	chr18	56863663	58638709	1775047	-13/2	BP:__TE_TE
inso	NA	760_0/754_0	6/6	2/6	chr1	188536417	188536417	-1	chr18	10930057	10930057	1	-20/-20	BP:TE_TE__
inso	NA	682_0/589_0	7/9	2/9	chr1	246112340	246112340	-1	chr18	66074916	66074916	1	-2/-2	BP:__TE_TE
inso	NA	504_0/497_0	15/16	1/10	chr11	80370988	80370988	-1	chr18	66199314	66199314	1	-34/-34	BP:TE_TE__
inso	NA	831_0/752_0	5/6	5/6	chr11	89609742	89609742	-1	chr18	20675062	20675062	1	-28/-28	BP:__TE_TE
inso	NA	554_0/544_0	12/12	2/5	chr13	89361761	89361761	-1	chr18	40628245	40628245	1	-10/-10	BP:TE_TE__
inso	NA	482_0/477_0	18/18	2/10	chr13	94154988	94154988	-1	chr18	31841819	31841819	1	-20/-20	BP:__TE_TE
inso	NA	552_0/543_0	12/12	2/10	chr13	94154988	94154988	-1	chr18	33366658	33366658	1	-24/-24	BP:___
inso	NA	572_0/582_0	11/10	2/10	chr13	94154988	94154988	-1	chr18	33414001	33414001	1	-25/-25	BP:__TE_TE
inso	NA	832_0/823_0	5/5	2/10	chr13	94154988	94154988	-1	chr18	33840728	33840728	1	-31/-31	BP:___
inso	NA	547_0/533_0	12/12	2/6	chr14	88453875	88453875	-1	chr18	47355345	47355345	1	-4/-4	BP:TE_TE_TE_TE
inso	NA	460_0/459_0	21/21	1/7	chr15	80290558	80290558	-1	chr18	10204472	10204472	1	-35/-35	BP:___
inso	NA	548_0/534_0	12/12	1/7	chr15	80290577	80290577	-1	chr18	10373958	10373958	1	-39/-39	BP:__TE_TE
inso	NA	549_0/535_0	12/12	1/10	chr18	10727847	10727847	-1	chr8	122638945	122638945	1	-2/-2	BP:___
inso	NA	637_0/718_0	8/6	1/3	chr18	41924326	41924326	-1	chr3	52215635	52215635	1	3/3	BP:___
inso	NA	503_0/501_0	15/15	1/6	chr18	50103558	50103558	-1	chr20	43152025	43152025	1	-2/-2	BP:__TE_TE
inso	NA	553_0/565_0	12/11	2/18	chr18	68315631	68315631	-1	chr8	110628783	110628783	1	-12/-12	BP:__TE_TE
inso	NA	445_0/443_0	25/25	2/10	chr18	70013527	70013527	-1	chr3	89469319	89469319	1	-30/-30	BP:__TE_TE
inso	NA	465_0/463_0	20/20	2/10	chr18	70123134	70123134	-1	chr3	89469319	89469319	1	-27/-27	BP:TE_TE_TE_TE
inso	NA	517_0/481_0	14/18	2/1	chr2	12668402	12668402	-1	chr18	14194320	14194320	1	-25/-25	BP:TE_TE__
inso	NA	827_0/639_0	5/8	1/5	chr2	32227442	32227442	-1	chr18	50012628	50012628	1	-40/-40	BP:TE_TE_TE_TE
inso	NA	469_0/468_0	19/19	1/7	chr20	43152025	43152025	-1	chr18	50430831	50430831	1	-29/-29	BP:TE_TE_TE_TE
inso	NA	746_0/678_0	6/7	6/2	chr20	56114021	56114021	-1	chr18	26134579	26173402	38824	-9/-21	BP:___
inss	NA	514_0/518_0	14/14	1/4	chr1	28801685	28801685	-1	chr18	77358872	77358872	1	-31/-31	BP:TE_TE__
inss	NA	607_0/610_0	9/9	1/13	chr11	89609742	89609742	-1	chr18	27827051	27827051	1	-28/-28	BP:___
inss	NA	522_0/519_0	13/14	1/4	chr13	72703887	72703887	-1	chr18	77358872	77358872	1	-9/-9	BP:TE_TE__
inss	NA	680_0/609_0	7/9	2/8	chr13	72703887	72703887	-1	chr18	79263837	79263837	1	-30/-30	BP:TE_TE__
inss	NA	758_0/684_0	6/7	1/10	chr18	10930057	10930057	-1	chr2	222878521	222878521	1	-25/-25	BP:__SR_SR
inss	NA	689_0/633_0	7/8	1/4	chr18	26134579	26134579	-1	chr9	10858712	10858712	1	-20/-20	BP:__TE_TE
inss	NA	834_0/824_0	5/5	1/2	chr18	26134579	26134579	-1	chrX	145458414	145458414	1	-41/-41	BP:___
inss	NA	608_0/581_0	9/10	1/5	chr18	26173402	26173402	-1	chrX	145458414	145458414	1	-29/-29	BP:___
inss	NA	611_0/600_0	9/9	11/1	chr18	45311409	45311409	-1	chr19	3982132	3982132	1	-33/-33	BP:TE_TE__
inss	NA	635_0/790_0	8/5	1/6	chr18	50012628	50012628	-1	chr2	127641529	127641529	1	-6/-6	BP:TE_TE__
inss	NA	470_0/479_0	19/18	1/5	chr18	61462044	61462044	-1	chr7	54378346	54378346	1	-18/-18	BP:__TE_TE
inss	NA	829_0/717_0	5/6	1/5	chr18	61462044	61462044	-1	chr7	143940165	143940165	1	-2/-2	BP:___
inss	NA	546_0/532_0	12/12	1/6	chr18	62500026	62500026	-1	chr7	54378346	54378346	1	-6/-6	BP:SR_SR_TE_TE
inss	NA	527_0/508_0	13/14	1/4	chr18	68315630	68315630	-1	chr3	89469299	89469299	1	1/1	BP:__TE_TE
inss	NA	545_0/531_0	12/12	1/10	chr18	69680373	69680373	-1	chr3	89469300	89469300	1	-37/-37	BP:TE_TE_TE_TE
inss	NA	688_0/686_0	7/7	1/6	chr18	69714982	69714982	-1	chr3	89469319	89469319	1	-38/-38	BP:__TE_TE
inss	NA	499_0/492_0	16/16	1/16	chr18	71159390	71159390	-1	chr3	8842566	8842566	1	-4/-4	BP:__TE_TE
inss	NA	440_0/439_0	28/28	1/17	chr18	72319702	72319702	-1	chr3	8842566	8842566	1	-9/-9	BP:__TE_TE
inss	NA	505_0/502_0	15/15	1/17	chr18	72711991	72711991	-1	chr3	8842566	8842566	1	-21/-21	BP:TE_TE_TE_TE
inss	NA	449_0/453_0	24/23	1/5	chr20	25658177	25658177	-1	chr18	49300016	49300016	1	-6/-6	BP:TE_TE__
inss	NA	551_0/571_0	12/11	1/5	chr20	25658177	25658177	-1	chr18	49418901	49418901	1	-26/-26	BP:TE_TE_SR_SR
inss	NA	451_0/450_0	24/24	1/7	chr20	43152004	43152004	-1	chr18	54763252	54763252	1	-22/-22	BP:TE_TE_TE_TE
inss	NA	641_0/685_0	8/7	2/7	chr3	192474316	192474316	-1	chr18	10930055	10930055	1	-29/-29	BP:TE_TE__
transl_inter	NA	570	11	8	chr1	66753472	-1	chr18	70444674	-1	-50	BP:TE_
transl_inter	NA	791	5	5	chr1	107988319	1	chr18	35980668	-1	-42	BP:_TE
transl_inter	NA	755_0	6	1	chr10	64140812	1	chr18	10930055	-1	-29	BP:TE_
transl_inter	NA	523_0	13	2	chr10	80363269	-1	chr18	66130374	-1	3	BP:_TE
transl_inter	NA	512	14	15	chr11	25649769	-1	chr18	75846882	-1	-43	BP:TE_
transl_inter	NA	789_0	5	11	chr12	11624029	-1	chr18	56863664	-1	-8	BP:_TE
transl_inter	NA	513	14	10	chr13	18843690	-1	chr18	14216367	-1	-45	BP:TE_TE
transl_inter	NA	590	9	13	chr14	66621359	1	chr18	22506567	-1	-45	BP:_
transl_inter	NA	632	8	20	chr15	66008014	-1	chr18	79263837	-1	-15	BP:TE_
transl_inter	NA	536	12	13	chr18	7650940	-1	chr3	87425411	1	-45	BP:TE_
transl_inter	NA	833	5	6	chr18	10930057	1	chr9	36676099	-1	-13	BP:_TE
transl_inter	NA	757_0	6	1	chr18	26134579	1	chr20	56114021	1	-21	BP:_
transl_inter	NA	759	6	1	chr18	26134579	1	chr3	18033189	-1	-21	BP:_TE
transl_inter	NA	826	5	6	chr18	26134579	-1	chr4	27670343	-1	-11	BP:_
transl_inter	NA	716	6	7	chr18	39236217	-1	chr3	166976579	1	-44	BP:TE_
transl_inter	NA	725	6	20	chr18	39236217	-1	chr8	52954930	1	29	BP:TE_
transl_inter	NA	452	24	7	chr18	43989113	1	chr19	3982132	-1	-15	BP:_
transl_inter	NA	436	30	8	chr18	50862264	1	chr20	43152025	-1	-8	BP:TE_TE
transl_inter	NA	480_0	18	8	chr18	56863663	1	chr6	147769543	1	2	BP:TE_
transl_inter	NA	462	20	6	chr18	58742568	-1	chr7	143940183	-1	-6	BP:_
transl_inter	NA	657	7	19	chr18	60932220	-1	chr5	142136757	-1	-44	BP:_TE
transl_inter	NA	828_0	5	3	chr18	66074917	1	chr7	48942925	1	-39	BP:TE_
transl_inter	NA	604_0	9	2	chr18	70444674	-1	chr7	48942943	1	-33	BP:_
transl_inter	NA	605_0	9	6	chr18	70444674	-1	chr6	24593189	1	-43	BP:_
transl_inter	NA	458	21	12	chr18	75846882	-1	chr5	140515479	1	-43	BP:_
