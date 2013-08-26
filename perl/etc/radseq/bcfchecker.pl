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
	my $t = "$DP $GT,$GQ $PL,$SP";
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
	print "$CHROM,$POS $REF,$ALT $QUAL $INFO, DP GT,GQ,SP\n";
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

