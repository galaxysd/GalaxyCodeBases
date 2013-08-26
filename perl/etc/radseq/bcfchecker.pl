#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
Purpose: Read bcf, get tped for p-link
Notes: rad2marker is deprecated.
=cut
use strict;
use warnings;
use Data::Dump qw(ddx);
use Galaxy::IO;
use Galaxy::SeqTools;
use Data::Dumper;

die "Usage: $0 <indexed bcgv bcf> <chr> <position>\n" if @ARGV<3;
my $bcfs=shift;
my $chr=shift;
my $pos=shift;

open L,'<','/share/users/huxs/work/tiger2/bcf/tiger2.tfam' or die;
my (@Case,@Control,%Pheno);
while (<L>) {
	next if /^(#|$)/;
	#print OF $_;
	chomp;
	my ($family,$ind,$P,$M,$sex,$pho) = split /\t/;
	if ($pho == 1) {
		push @Case,$ind;
	} elsif ($pho == 2) {
		push @Control,$ind;
	} else { die; }	# $pho can only be 1 or 2
	$Pheno{$ind} = $pho;	# disease phenotype (1=unaff/ctl, 2=aff/case, 0=miss)
}

open IN,'-|',"bcftools view $bcfs $chr:$pos-$pos" or die "Error opening output $bcfs:$!\n";
my (@Samples,@Parents);
while (<IN>) {
	next if /^##/;
	chomp;
	my @data = split /\t/;
	if ($data[0] eq '#CHROM') {
		@Samples = map {
			s/_2$//;
			my $t = $_;
		} splice @data,9;
		@Parents = grep(/^GZXJ0(34|35|39)$/,@Samples);	# P P g g g P g g g(9) GGGGGGGGGG(10)
		last;
	}
}

sub FormatIt($) {
	my $str = $_[0];
	my ($GT,$PL,$DP,$SP,$GQ) = split /:/,$str;
	my $t = "$DP\t$GT:$GQ\t$PL:$SP";
	$t = "\t\t\t\t\t$t" if $DP == 0;
	return $t;
}

while (<IN>) {
	next if /^#/;
	chomp;
	my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, @data) = split /\t/;
	print "$_\n";
	my %Dat;
	@Dat{@Samples} = @data;
	#ddx \%Dat;
	print "$CHROM,$POS $REF,$ALT $QUAL $INFO, DP GT:GQ PL:SP\n";
	print "Parents\n";
	for (@Parents) {
		print FormatIt($Dat{$_}),"\t$_\n";
	}
	print "Control\n";
	for (@Control) {
		print FormatIt($Dat{$_}),"\t$_\n";
	}
	print "Case\n";
	for (@Case) {
		print FormatIt($Dat{$_}),"\t$_\n";
	}	
}
close IN;

__END__
perl bcfchecker.pl snow_white_000000_210210.bcgv1.bcf scaffold97 5273427

scaffold97 1699193,5985068
scaffold1457 1416415,6191481

CORIN
scaffold1457	2643215,2892269

grep \, genes1457.candidateDNA.fa
cat genes1457.log

1457:
156935 x

345865 x
345868 x
345876 ?	SCFD2
345882 ?
345887 ?
345890 ?

734940 ?
1323104 ?
1332800 ?
1604103 ?
1955320 ?
2375977 ?
2379975 x
2434696 x
2477112 ?
2488325 x
2550691 ?
2568107 Y
2635027 ?

2725459 ?
2725471 ?
2798540 x	W317R	CORIN
2829955 ?	H587Y
2866888 x

2896448 ?
2946009 ?
3485780 ?
3695806 ??
5396351 ?

97:
1151252 x
1189163 ?
1189168 ?
1429506 x
1798413 ?
1940991 ?
2163261 ?
2165439 x
2191209 ?
2191308 ?
2299928 x
2309017 ?
2354904 ?
2380816 ??
2569841 x
2584490 Y?
2590189 ?
2741630 ?
2785256 ?
3347491 ?
3347621 ?
3472317 ?
3472380 ?
3552459 ?
3556171 ?
3556398 ?
3612742 ?
3617673 ?
3622659 x
3925885 x
3929632 ??
3929668 ?
3990107 x
3993753 ?
4078789 ?
4300974 x
4340409 ?
4786540 ?
4799288 ?
4799292 ?
4969678 ?
5474455 ?
5775713 x
