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

die "Usage: $0 <bcgv bcf> <min Sample Count> <out>\n" if @ARGV<3;
my $bcfs=shift;
my $minSampleCnt=shift;	# 16 for 16 samples in paper
my $outfs=shift;

my (%Stat,$t);
open OV,'|-',"gzip -9c > $outfs.vcf.gz" or die "Error opening $outfs.vcf.gz: $!\n";

my (%Pheno,@tfamSamples,%tfamDat,%tfamSampleFlag,%inFamily,%ISinFamily);

my $th = openpipe('bcftools view --types snps -m2',$bcfs);	# -I	skip indels, para changed in v1.1; '-m2 -M2' for biallelic sites.
warn "[!]Begin.\n";
my @Samples;
while (<$th>) {
	print OV $_;
	next if /^##/;
	chomp;
	my @data = split /\t/;
	if ($data[0] eq '#CHROM') {
		@Samples = map {
			if (/^(\w+)_all\.bam$/) {
				$_ = $1;
			} elsif (/\.\/bam0\//) {
				my $t=(split /_/)[-1]; $t=(split /\./,$t)[0]; $_=$t;
			} else {
				#my $t=(split /\//)[-1];$t=~s/_cut//g;$t=~s/-/./g; $_=join('_',(split /\./,$t)[-5,-6]);
				s/_2$//;
				my $t = $_;
			}
		} splice @data,9;
		last;
	}
}
while (<$th>) {
	if (/^#/) {
		print OV $_;
		next;
	}
	chomp;
	my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, @data) = split /\t/;
	++$Stat{'VCF_In'};	# also as rs#
	my @Bases = split /\s*,\s*/,"$REF, $ALT";
	my @groups = split(/\s*;\s*/, $INFO);
	my (%INFO,$name);
	for my $group (@groups) {
		my ($tag,$value) = split /=/,$group;
#warn "- $group -> $tag,$value\n";
		my @values = split /,/,$value;
		if (@values == 1) {
			$INFO{$tag}=$values[0];
		} else {
			$INFO{$tag}=\@values;
		}
	}
	my (%GT,%GTcnt);
	my @FMT = split /:/,$FORMAT;
	for my $s (@Samples) {
		my $dat = shift @data or die "bcf file error.";
		my @dat = split /:/,$dat;
		for my $i (@FMT) {
			$GT{$s}{$i} = shift @dat;
		}
	}
	my $SPcnt = 0;
	my $GQok = 1;
	my (%GTitemCnt,$Mut,@plinkGT);
	for (@Samples) {
		if ($GT{$_}{'GQ'} < 20) {
			$GQok = 0;
			last;
		}
		if ($GT{$_}{'DP'} > 0) {
			my $gt = $GT{$_}{'GT'};
			++$GTcnt{$gt};
			my @GT = split /[\/|]/,$gt;
			++$SPcnt;
			$GT{$_}{'GTp'} = join ' ',map($Bases[$_], @GT);
			#$GT{$_}{'O_K'} = 1;
		} else {
			$GT{$_}{'GTp'} = '0 0';
			#$GT{$_}{'O_K'} = 0;
		}
		push @plinkGT,$GT{$_}{'GTp'};
	}
	#$Mut = $Bases[$Mut];
=pod
++$Stat{'GTcnt'}{$INFO{'FQ'} <=> 0}{scalar(keys %GTcnt)};
ddx $Stat{'GTcnt'};
ddx $CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO,\%INFO,\%GT if scalar(keys %GTcnt) > 1 and $INFO{'FQ'} < 0 and $SPcnt>2;
# rad2marker.pl:135: {
#   -1 => { 1 => 2850, 2 => 526, 3 => 37 },
#   1  => { 1 => 8, 2 => 2507, 3 => 1792 },
# }
=cut
#warn "$CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO\n";
#ddx $CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO,\%GTcnt,\%INFO,\%GT,\%GTitemCnt,$Mut;
	if ( $GQok == 0 or $SPcnt < $minSampleCnt) {	# No 'FQ' in 'INFO' for v1.1
		++$Stat{'VCF_Skipped'};
		#warn $_,"\n";
		next;
	}
	++$Stat{'Marker_Out'};
	print OV $_,"\n";
}
close $th;
close OV;

ddx \%Stat;

print "[Prepare $outfs.phe.] And then:\np-link --tfile $outfs --reference-allele $outfs.MinorAllele --fisher --out ${outfs}P --model --cell 0 --allow-no-sex\n--pheno $outfs.phe --all-pheno [screen log will be saved by p-link itselt]\n";
__END__

gatk -T CoveredByNSamplesSites -R ../ref/Felis_catus80_chr.fa -minCov 1 -percentage 0.5 -V filtered2.vcf -out ttt.out -l DEBUG -et STDOUT -K /opt/jar/huxs_pku.edu.cn.key
bcftools query -f '%CHROM,%POS\t%REF|%ALT|%QUAL\t%FILTER|%DP[\t%SAMPLE=%GT,%DP,%GQ]\n' -T <(sed 's/:/\t/' ttt.out) filtered2.vcf |les

bcftools view -s^FCAP114 -m2 mpileup_20150321HKT071334.vcf.gz \
 | bcftools norm -Df ../ref/Felis_catus80_chr.fa -c e -m+both \
 | bcftools filter -sLowQual -e'%QUAL<10' \
 | bcftools filter -m+ -sDepthHigh -e'DP>650' \
 | bcftools filter -m+ -sDepthLow -e'DP<2' \
 | bcftools filter -m+ -sBadSites -e'%QUAL<10 && RPB<0.1' \
 | tee >(bcftools view -Oz -o filter_20150321HKT071334.vcf.gz) \
 | bcftools view -f .,PASS -Oz -o filtered.vcf.gz &

bcftools query -f '%CHROM,%POS\t%REF|%ALT|%QUAL\t%DP[\t%SAMPLE=%GT,%DP]\n' filtered.vcf.gz |les

./bcf2ped.pl ~/work/catbtail/merged/kinkcats.tfam ~/work/catbtail/merged/mpileup_20150321HKT071334.vcf.gz 18 D xxxxx 2>xxxxx
