#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);
use Galaxy::IO;
#use Galaxy::SeqTools;
#use Data::Dumper;

my $WinSize=10,000;
my $RegionBegin=151586958;
my $RegionEnd=152939134;
my $Regions='gi|753572091|ref|NC_018727.2|:'."$RegionBegin-$RegionEnd";

die "Usage: $0 <mpileup bcf> <tfam file> <out>\n" if @ARGV<3;
my $bcfs=shift;
my $tfamfs=shift;
my $outfs=shift;

my (@tfamSamples,%tfamSamplePheno,%inFamily,@CaseS,@ControlS);
open F,'<',$tfamfs or die $!;
while (<F>) {
	next if /^(#|$)/;
	chomp;
	my ($family,$ind,$P,$M,$sex,$pho) = split /\t/;
	if ($pho == 1 or $pho == 2) {
		$tfamSamplePheno{$ind} = $pho;	# disease phenotype (1=unaff/ctl, 2=aff/case, 0=miss)
		if ($pho == 2) {
			push @CaseS,$ind;
		} else {
			push @ControlS,$ind;
		}
	} else { die; }	# $pho can only be 1 or 2
	push @tfamSamples,$ind;
	push @{$inFamily{$family}},$ind;
}
close F;

my $cmd = 'bcftools query -f \'%CHROM\t%POS\t%REF\t%DP[\t%SAMPLE=%DP]\n\' -r \'' . "$Regions' -s " . join(',',@CaseS,@ControlS);
warn "$cmd $bcfs |less -S\n";
my $fh = openpipe($cmd,$bcfs);
while (<$fh>) {
	chomp;
	my ($chr,$pos,$ref,$adp,@sampleDat) = split /\t/;
	#print "$chr,$pos,$ref:@sampleDat\n";
	for (@sampleDat) {
		my ($id,$sdp) = split /=/,$_;
	}
}
close $fh;

__END__
./bcfdeplot.pl mpileup_20150402HKT165931.bcf outA13.tfam outA13.depth10k



bcftools query -f '%CHROM,%POS\t%DP[\t%SAMPLE=%DP]\n' mpileup_20150402HKT165931.bcf|les

http://stackoverflow.com/questions/8714355/bash-turning-multi-line-string-into-single-comma-separated
grep -P '\t2$' outA13.tfam|awk -vORS=, '{print $2}'| sed 's/,$/\n/'
grep -P '\t2$' outA13.tfam|awk '{print $2}'| paste -d, -s

bcftools query -f '%CHROM,%POS\t%DP[\t%SAMPLE=%DP]\n' -s `grep -P '\t2$' outA13.tfam|awk '{print $2}'| paste -d, -s` mpileup_20150402HKT165931.bcf|les
