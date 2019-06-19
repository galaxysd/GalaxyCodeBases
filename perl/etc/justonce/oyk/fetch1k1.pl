#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 1.0.0 @ 20180516
=cut
use strict;
use warnings;

use Data::Dump qw(ddx);

my (%Pos19d,%Pos38d);
open IN,'<','nippt5427.tsv' or die $!;
while (<IN>) {
	next if /^#/;
	chomp;
	my ($rs,$chr19,$pos19,$chr38,$pos38,@d) = split /\t/,$_;
	$chr19 =~ s/^chr//;
	$Pos19d{$chr19}{$pos19} = $rs;
	$Pos38d{$chr38}{$pos38} = $rs if $chr38 ne 'unlocalized';
}
#ddx \%Pos19d;

my $hg19fa = '/share/FGI2017B/pub/Ref/Human/hs37d5.fna.bgz';
my $hg38fa = '/share/FGI2017B/pub/Ref/Human/GRCh38_no_alt_analysis_set.fna.bgz';
my @ChrIDs19 = (1..22,'X','Y','MT');
my @ChrIDs38 = map {"chr$_"} (1..22,'X','Y','M');

sub dealfa($$$$) {
	my ($PosH,$ChrIDsA,$faf,$outf) = @_;
	#ddx ($outf,$faf); ddx $ChrIDsA;
	#ddx $PosH;
	open O,'>',$outf or die $?;
	for my $chr (@$ChrIDsA) {
		my @Poses = sort { $a <=> $b } keys %{$PosH->{$chr}};
		for my $p (@Poses) {
			my $left = $p - 500;
			my $right = $p + 500;
			die if $left < 1;
			open( SEQ,"-|","samtools faidx $faf $chr:$left-$right") or die $!;
			my $head = <SEQ>;
			if ($head ne ">$chr:$left-$right\n") {
				die;
			}
			my $theSeq;
			while (<SEQ>) {
				chomp;
				s/\s//g;
				$theSeq .= $_;
			}
			my $top = substr $theSeq,0,500;
			my $base = substr $theSeq,500,1;
			my $bottom = substr $theSeq,501,500;
			print O ">$PosH->{$chr}{$p} $chr:$p $left-$right\n$top\n$base\n$bottom\n";
		}
	}
	close O;
}
dealfa(\%Pos19d,\@ChrIDs19,$hg19fa,'nippt5427.hg19.fa');
dealfa(\%Pos38d,\@ChrIDs38,$hg38fa,'nippt5427.GRCh38.fa');


__END__
$ grep -v GL hs37d5.fna.bgz.fai
1	249250621	52	60	61
2	243199373	253404903	60	61
3	198022430	500657651	60	61
4	191154276	701980507	60	61
5	180915260	896320740	60	61
6	171115067	1080251307	60	61
7	159138663	1254218344	60	61
8	146364022	1416009371	60	61
9	141213431	1564812846	60	61
10	135534747	1708379889	60	61
11	135006516	1846173603	60	61
12	133851895	1983430282	60	61
13	115169878	2119513096	60	61
14	107349540	2236602526	60	61
15	102531392	2345741279	60	61
16	90354753	2449981581	60	61
17	81195210	2541842300	60	61
18	78077248	2624390817	60	61
19	59128983	2703769406	60	61
20	63025520	2763883926	60	61
21	48129895	2827959925	60	61
22	51304566	2876892038	60	61
X	155270560	2929051733	60	61
Y	59373566	3086910193	60	61
MT	16569	3147273397	70	71
NC_007605	171823	3153506529	60	61
hs37d5	35477943	3153681224	60	61

$ egrep -v 'chrUn|_random' GRCh38_no_alt_analysis_set.fna.bgz.fai
chr1	248956422	112	70	71
chr2	242193529	252513167	70	71
chr3	198295559	498166716	70	71
chr4	190214555	699295181	70	71
chr5	181538259	892227221	70	71
chr6	170805979	1076358996	70	71
chr7	159345973	1249605173	70	71
chr8	145138636	1411227630	70	71
chr9	138394717	1558439788	70	71
chr10	133797422	1698811686	70	71
chr11	135086622	1834520613	70	71
chr12	133275309	1971537157	70	71
chr13	114364328	2106716512	70	71
chr14	107043718	2222714743	70	71
chr15	101991189	2331287770	70	71
chr16	90338345	2434736088	70	71
chr17	83257441	2526365093	70	71
chr18	80373285	2610812039	70	71
chr19	58617616	2692333639	70	71
chr20	64444167	2751788762	70	71
chr21	46709983	2817153685	70	71
chr22	50818468	2864531079	70	71
chrX	156040895	2916075638	70	71
chrY	57227415	3074345836	70	71
chrM	16569	3132390908	70	71
chrEBV	171823	3144056708	70	71
