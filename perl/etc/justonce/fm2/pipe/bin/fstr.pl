#!/usr/bin/env perl

use strict;
use warnings;

#use Data::Dump qw(ddx);

my $s1HomoLevel = 0.8;
my $s2HomoLevela = 0.75;
my $s2HeteLevela = 0.3;
my $s2HeteLevelb = 0.2;
my $minAllelDepth = 100;

my @STRids = qw[
CSF1PO D10S1248 D10S1435 D12S391 D13S317 D15S659 D16S539 D19S253 D19S433 D1S1656 D21S11 D22S1045 D2S1338 D2S441 D3S1358 D3S3045 D5S818 D6S1043 D6S477 D7S820 D8S1132 D8S1179 DXS6795 FGA PentaE SE33 TH01 TPOX Y_GATA_A10 Y_GATA_H4 vWA DYF387 DYF404S1 DYS19 DYS385 DYS388 DYS390 DYS391 DYS392 DYS393 DYS437 DYS438 DYS439 DYS444 DYS447 DYS448 DYS449 DYS456 DYS458 DYS459 DYS460 DYS481 DYS510 DYS518 DYS522 DYS527 DYS531 DYS533 DYS549 DYS557 DYS570 DYS576 DYS593 DYS596 DYS627 DYS635 DYS645
];
my %STRid;
$STRid{$_}++ for (@STRids);

my (%StrZones,%StrStat);
<DATA>;
while(<DATA>) {
	chomp;
	my ($id,undef,undef,$loc,$length) = split /\t/;
	$StrZones{$id} = [$loc,$length];
	#push @{$StrStat{$length}},$id;
}
#ddx \%StrZones;
#ddx \%StrStat;

my %StrDat;
while(<>) {
	chomp;
	my ($aID,$depth) = split /\s+/;
	my ($strID,$strAllele) = split /\-/,$aID;
	#print join("\t",$strID,$strAllele,$depth),"\n";
	$StrDat{$strID}{$strAllele} = $depth;
}
#ddx \%StrDat;

sub SumArryItr { my $agg = 0; $agg += $_ for @_;  return $agg }

for my $k (sort { my $x=($a=~/^DY/);my $y=($b=~/^DY/);$x-$y ||$a cmp $b} keys %STRid) {
	my $strLen = 0;
	$strLen = $StrZones{$k}->[1] if exists $StrZones{$k};
	unless (exists $StrDat{$k}) {
		print "STR\t$k: Unknown\n";
		push @{$StrStat{$strLen}},'N';
	} else {
		my %oneSTR = %{$StrDat{$k}};
		#ddx \%oneSTR;
		my @Allels = sort {$oneSTR{$b} <=> $oneSTR{$a}} keys %oneSTR;
		#ddx \@Allels;
		my $Sum = SumArryItr(values %oneSTR);
		#print $Sum,"\n---\n";
		my @GT;
		if ($oneSTR{$Allels[0]} <= 100) {
			@GT = ('NA');
			push @{$StrStat{$strLen}},'N';
		} elsif ($oneSTR{$Allels[0]} >= $Sum * $s1HomoLevel ) {
			@GT = ($Allels[0]);
			push @{$StrStat{$strLen}},$StrDat{$k}{$Allels[0]};
		} elsif ( ($oneSTR{$Allels[0]} >= $Sum * $s2HomoLevela) and ($oneSTR{$Allels[1]} < $Sum * $s2HeteLevelb) ) {
			@GT = ($Allels[0]);
			push @{$StrStat{$strLen}},$StrDat{$k}{$Allels[0]};
		} elsif ( ($oneSTR{$Allels[0]} >= $Sum * $s2HeteLevela) and ($oneSTR{$Allels[1]} >= $Sum * $s2HeteLevelb) ) {
			@GT = ($Allels[0],$Allels[1]);
			push @{$StrStat{$strLen}},($StrDat{$k}{$Allels[0]}+$StrDat{$k}{$Allels[1]});
		} else {
			@GT = ('UnknownE');
			push @{$StrStat{$strLen}},'N';
		}
		#ddx \@GT;
		print "STR\t$k: ",join(",",@GT),"\n";
	}
}

#ddx \%StrStat;
for my $plen (sort {$a <=> $b} keys %StrStat) {
	next if $plen == 0;
	my @deps = @{$StrStat{$plen}};
	my ($cnt,$sum,$avg)=(0,0,'NA');
	for (@deps) {
		next if $_ eq 'N';
		$sum += $_;
		$cnt += 1;
	}
	if ($cnt) {
		$avg = sprintf("%.1f", $sum/$cnt);
	}
	print join("\t",'DEPTH',$plen,$avg,@deps),"\n";
}

__DATA__
ID	F	R	Location	Product_Size
AMEL	TTGGCTCACGTGACTCCAACCAGA	ACCCCTTTGAAGTGGTACCAGAGCA	chrY:6868493-6868770	278
CSF1PO	CCTGGTGCACACTTGGACAGCATT	TACAGAGGAGGCACTTCGTGGTGG	chr5:150076268-150076544	277
D10S1248	AGCAAACCTGAGCATTAGCCCCAG	TAAGTGCAGTGCTTGGCAAAGAGC	chr10:129294172-129294435	264
D10S1435	TGCTAGCACGTTGGGTTTCCTGAC	AATCAGCTGCTTGGTGGTGTGCAC	chr10:2201068-2201372	305
D12S391	TCAACAGGATCAATGGATGCATAGG	ACTCCAGGTTCTCAGGCCTTCCAC	chr12:12296994-12297283	290
D13S317	TTGGGATGGGTTGCTGGACATGGTA	AACTTGGGTTGAGCCATAGGCAGC	chr13:82147897-82148130	234
D15S659	CGTCTTCCCAACATAACATATTGCTT	AGGTGGGAGAGACAGACGGATACCT	chr15:46081870-46082178	309
D16S539	TGCCAGATGCTCGTTGTGCACAAA	GCGTTTGTGTGTGCATCTGTAAGC	chr16:86352550-86352807	258
D18S51	GAGTTCTTGAGCCCAGAAGGTTA	TGGTGTGTGGAGATGTCTTACAAT	chr18:63281578-63281811	234
D19S253	ACCTCTAAAGCCTGCAGTGGGACC	TTTGGACCCAGACGGAATTGCACC	chr19:15617389-15617698	310
D19S433	GGTGCACCCATTACCCGAATAAAA	GGTTGAGGCTGCAAAAAGCTATA	chr19:29926188-29926368	181
D1S1656	TTCCACCGCAGCACAAAACTCGTT	TTCAAGCCTGTGTTGCTCAAGGGT	chr1:230769446-230769720	275
D21S11	GGCTTCAGACTTGGACAGCCACAC	GTCAATGTTCTCCAGAGACAGACT	chr21:19181852-19182156	305
D22S1045	CGTTGGAATTCCCCAAACTGGCCA	TAGTGACAGTGCCTGTGCCCAAGT	chr22:37140181-37140461	281
D2S1338	GAGGCCCTTGTCAGTGTTCATGCC	ATAATCCAGCTGTGGGAGGGAGCC	chr2:218014756-218015005	250
D2S441	ACCCACGGCCAGAAAGTTGGGTAA	TTGGAGCTAAGTGGCTGTGGTGTT	chr2:68011767-68012024	258
D3S1358	AGCTATTCCCAGGTGAGGACTGCA	TTCCCCCACTGCAGTCCAATCTGG	chr3:45540569-45540850	282
D3S3045	TGGCATAAATCAGGACCTCTCTGGGC	ACAGCACAGATCCAATGTCAGCAA	chr3:107271085-107271381	297
D5S818	ACCTCTCCCATCTGGATAGTGGACC	TGTGACAAGGGTGATTTTCCTCTTTGG	chr5:123775371-123775647	277
D6S1043	AGCAGCCTACCATGTTTTGAAGGC	AGTGTGCAAGGATGGGTGGATCAA	chr6:91740036-91740310	275
D6S477	TGTCACAGGGCTGATGAGGTGAAAT	TGTGTCTCAGGTAGCAGCAGGACT	chr6:6140356-6140665	310
D7S820	TCCTCATTGACAGAATTGCACCAAA	TGTTGGTCAGGCTGACTATGGAGT	chr7:84160138-84160366	229
D8S1132	TCCTTGTTTCCTCATTTTATTTCGGTCCC	CTCTCTCCCTCTCTCTTTCGAGCA	chr8:106316564-106316799	236
D8S1179	GATCCTTGGGGTGTCGCTTTTCTGG	CGTATCCCATTGCGTGAATATGCCT	chr8:124894720-124894984	265
DXS6795	TGCATCCATCCCCTAAACCTCTCA	TACCGGTGGGATTCTGGGTGTCTG	chrX:23226378-23226684	307
DYF387S1a	GAGCTAGATTCCATTTTACCCCTAACA	CTGGTGCCACAGTGTGAGAAGTGT	chrY:23785305-23785554	250
DYF387S1b	AGTGTGAGAAGTGCTACCACAGTT	GGGTGACAGAGCTAGATTCCATTTTACCC	chrY:25884546-25884788	243
DYF404S1a	CAAAGGGCTTAAGAAATTTCAACGC	CCCAGGATTGAGAGGCTGCAGACT	chrY:23807903-23808136	234
DYF404S1b	CAAAGGAGCCCAGGATTGAGAGGC	TCAAAGGGCTTAAGAAATTTCAACGCA	chrY:25861945-25862179	235
DYS19	GCCATGGCCATGTAGTGAGGACAA	TGACAAGCCCAAAGTTCTTAACATTC	chrY:9684322-9684621	300
DYS385a	TGGTAAGGGCTGCTGACCAGATTTCT	AGAGCTAGACACCATGCCAAACAACAA	chrY:18639588-18639926	339
DYS385b	GACACCATGCCAAACAACAACAAAGAAAAGA	GGGCTGCTGACCAGATTTCTTTCTGA	chrY:18680469-18680806	338
DYS388	AGCTTTTAGTGAGCCAGATCGCACC	GGTCTGGGCAGCGAGTTTAAGGGA	chrY:12635538-12635804	267
DYS389I	AGATGCGAGAACACCCATCCCAGT	CCTACTTCTGTATCCAACTCTCATCTGT	chrY:12500293-12500535	243
DYS389II	GGAACACAATTATCCCTGAGTAGCA	ACCTACTTCTGTATCCAACTCTCATCT	chrY:12500361-12500656	296
DYS390	AGTCCTGAGACAGTGTATCCGCCAT	TTTGGGCCCTGCATTTTGGTACCC	chrY:15162942-15163204	263
DYS391	GCACAAGACACCCCACCACAGATT	TCCCTGGTTGCAAGCAATTGCCA	chrY:11981878-11982187	310
DYS392	GGGATCATTAAACCTACCAATCCCA	CCATCCATGTTGCTCCAAAGGACCC	chrY:20471951-20472227	277
DYS393	TCCTAATGTGGTCTTCTACTTGTGTC	GCCAGATAACGTGTGTGGAAGTGA	chrY:3263080-3263344	265
DYS437	TAGCTGGGACTATGGGCGTGAGTG	TGGGCATGCACTCACACCCATAGT	chrY:12346230-12346554	325
DYS438	GGGGAATAGTTGAACGGTAAACAGT	GAAAATATACCGTTTACTGGCCTGG	chrY:12825861-12826255	395
DYS439	ACTTCCTAGGTTTTCTTCTCGAGTTG	TGCCTGGCTTGGAATTCTTTTACCCA	chrY:12403375-12403612	238
DYS444	TGAAAGGTGTGAACCATTTGGCATGT	TGAGCCCATGCCATTCAAACTCACG	chrY:17114136-17114410	275
DYS447	TACCTGTGACCTGCAAGCTCCCAA	GGTGGTCACAGCATGGCTTGGTTT	chrY:13166712-13167004	293
DYS448	AGATCGCGAGACAGAAAGGGAGAT	TCCTCATATTTCTGGCCGGTCTGGA	chrY:22218885-22219128	244
DYS449	CTTGGAGTCTCTCAAGCCTGTTCT	GGTTGGACAACAAGAGTAAGACAGAA	chrY:8349865-8350173	309
DYS456	AGTTTTGGAACTCGGACTGGCTCAT	ACCCCATCAACTCAGCCCAAAACT	chrY:4402844-4403051	208
DYS458	GCAACAGGAATGAAACTCCAATGAAAG	GATTTCCTGACCTTGTGATCCAGC	chrY:7999799-7999983	185
DYS459a	ACCTTCCACCGTAAATAAACCCACCA	CACACCCCACTGCAGGCCAAAAAT	chrY:23932515-23932770	256
DYS459b	CCCCACTGCAGGCCAAAAATAGGT	ACCTTCCACCGTAAATAAACCCACCA	chrY:25737350-25737597	248
DYS460	TGTCCAGTAGTGATGCTGTGTCAC	GGAAAGTCAAGACAGTAGCAAGCACA	chrY:18888759-18889051	293
DYS481	TTTCTGCTAGGTGGGGCTTCCTGT	CCCACAACCCAAGAAGAGCCACAC	chrY:8558150-8558459	310
DYS510	TAGATCACAGGCCCTGTCGCAGAA	ATCCATCCATCCAACTGTTCTCCC	chrY:15187889-15188113	225
DYS518	GGGCAACACAAGTGAAACTGCTTCTCG	GCACTCTGGGCTCCAAGAATAGTGGA	chrY:15207958-15208284	327
DYS522	TGTGCATGGGGAAGCTCTGATTGA	AACCAGTGAGAGCCGGAACCTCAT	chrY:7547385-7547661	277
DYS527a	TCGCAAACATAGCACTTCAGCCCAG	AGCCACAACATAAGTAAGGTAGTTTTCT	chrY:23739457-23739803	347
DYS527b	AGCCACAACATAAGTAAGGTAGTTTTCT	ATCGCAAACATAGCACTTCAGCCCA	chrY:25930281-25930620	340
DYS531	TGTCTTTGTGGCTTTGCTTGGGCT	AGATACATTCAAGGGCCCTGCGGA	chrY:8598034-8598334	301
DYS533	TCTTCTACCTATCATCTTTCTAGCTAGCT	GCATCCAAAAGCTTTGAGTGTGGA	chrY:16281314-16281552	239
DYS549	AGTGTAAGCCAAACCCAAATATAGCA	TGGTGGCATAAGTGGTAATGTCCCC	chrY:19358153-19358450	298
DYS557	ACGCATCTGCCCAGCATGTGTTTT	ATGCACCTTGAGGGATGCCAAAGC	chrY:21072692-21072914	223
DYS570	CCTGGCTGTGTCCTCCAAGTTCCT	AGGGCTTCTAAGGGATGCAAGGTGT	chrY:6993166-6993367	202
DYS576	CCAAGCACAGTGGGTCATGCCAG	TCCTGGAGATGAAGGAGGAGATGGG	chrY:7185199-7185439	241
DYS593	GCAAGGGGCACATGCATGTCATCC	TGCTGTGTCCTACCTTTGAGAGCCT	chrY:16473740-16474046	307
DYS596	ATGTACAGGTCCAAAGGCAGCAGC	GGAAAATTTTGACAAGCCCAAAGTTC	chrY:8519474-8519772	299
DYS627	TAGGTGACAGCGCAGGATTCCATCT	GTTTCTGTGAGTCCACTGGAGACCT	chrY:8781890-8782284	395
DYS635	TCAGCTTGAGTGATGGACCAAGGC	TGCCCAATGGAATGCTCTCTTGGC	chrY:12258700-12258995	296
DYS643	AGCCATGCCTGGTTAAACTACTGTGC	AACTCAAGGGAAGGGTGGCAGGAT	chrY:15314106-15314368	263
DYS645	GGTTTTGGTTACGGGTGGCAATCA	GCACCTCACACTGAAGGATTGGCC	chrY:21165925-21166176	252
FGA	ATCCTCTGACACTCGGTTGTAGGT	TGGAAGGCTGCAGGGCATAACATT	chr4:154587623-154587930	308
Penta_D	AGGCATGGTGAGGCTGAAGTAGGA	TGCCTAACCTATGGTCATAACGATTTT	chr21:43636078-43636350	273
Penta_E	TTTGGGTTATTAATTGAGAAAACTCCTTAC	AAAATACATTTTACCAACATGAAAGGGTAC	chr15:96830980-96831369	390
SE33	TCTCCCCTACCGCTATAGTAACTTGC	GGTGCACGTCTGTAATTCCAGCTCC	chr6:88277104-88277491	388
TH01	TAAATGTGCCAGGGAGCCCAAGGT	CAAAGGGTATCTGGGCTCTGGGGT	chr11:2170931-2171173	243
TPOX	TCAGGGCTGTGATCACTAGCACCC	TTGACTCTACTGTCCTGGGCGCTC	chr2:1489583-1489875	293
Y_GATA_A10	ACCTACCTACCTATCCACCTGCCA	AGAGAGAGAGAAAGATAGAGATGGAAGGATAG	chrY:16606916-16607153	238
Y_GATA_H4	ACATAGCCCACTTGTTAAACAACTTAAC	GGGTTCTGAAGAGCTAAACAGAGACC	chrY:16631584-16631891	308
Yindel	ACCCAAATCAACTCAACTCCAGTGA	CTCTTGCAGCATTTTCAGTTAGCCT	chrY:13396761-13396981	221
vWA	AGGACAGATGATAAATACATAGGATGGATG	TAGTCTCCTACAATGTGCCGGGCA	chr12:5983937-5984227	291
