#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Data::Dump qw(ddx);
use Galaxy::IO;
use Galaxy::SeqTools;

die "Usage: $0 <tfam file> <bcgv bcf> <out>\n" if @ARGV<3;
my $tfamfs=shift;
my $bcfs=shift;
my $outfs=shift;

my (%Stat,$t);
open OP,'>',$outfs.'.tped' or die "Error opening $outfs.tped : $!\n";
open OM,'>',$outfs.'.MinorAllele' or die "Error opening $outfs.MinorAllele : $!\n";
open OD,'>',$outfs.'.dict' or die "Error opening $outfs.dict : $!\n";
if ($tfamfs ne $outfs.'.tfam') {
	open OF,'>',$outfs.'.tfam' or die "Error opening $outfs.tfam : $!\n";
}
open O,'>',$outfs.'.bcf2pedlog' or die "Error opening $outfs.bcf2pedlog : $!\n";
$t = "# In: [$bcfs], Out: [$outfs]\n";
print O $t;
print $t;

my (%Pheno,@tfamSamples);
open L,'<',$tfamfs or die;
while (<L>) {
	next if /^(#|$)/;
	print OF $_;
	chomp;
	my ($family,$ind,$P,$M,$sex,$pho) = split /\t/;
	$Pheno{$ind} = $pho;	# disease phenotype (1=unaff, 2=aff, 0=miss)
	push @tfamSamples,$ind;
}
close L;
close OF;

my $th = openpipe('bcftools view -I',$bcfs);	# -I	skip indels
my (@Samples,@Parents);
while (<$th>) {
	next if /^##/;
	chomp;
	my @data = split /\t/;
	if ($data[0] eq '#CHROM') {
		@Samples = map {
			if (/^(\w+)_all\.bam$/) {
				$_ = $1;
			} else {
				my $t=(split /\//)[-1];$t=~s/_cut//g;$t=~s/-/./g; $_=join('_',(split /\./,$t)[-5,-6]);
			}
		} splice @data,9;
		# ../5.bam_0000210210_merged/d1_4_merged.D4.JHH001.XTU.sort.rmdup.bam
		@Parents = grep(!/^GZXJ/,@Samples);
		last;
	}
}
print O "# Samples: [",join('],[',@Samples),"]\n# Parents: [",join('],[',@Parents),"]\n";
warn "Samples:\n[",join("]\n[",@Samples),"]\nParents: [",join('],[',@Parents),"]\n";
die "Samples in tfam and bcf not match !\n" if @tfamSamples != @Samples;
for (my $i = 0; $i < @Samples; $i++) {
	die "Samples in tfam and bcf not match !\n" if $tfamSamples[$i] ne $Samples[$i];
}	# http://stackoverflow.com/questions/2591747/how-can-i-compare-arrays-in-perl

while (<$th>) {
	next if /^#/;
	chomp;
	my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, @data) = split /\t/;
	++$Stat{'VCF_In'};	# also as rs#
	my @Bases = split /\s*,\s*/,"$REF, $ALT";
	my @groups = split(/\s*;\s*/, $INFO);
	my (%INFO,$name);
#	if ($groups[0] eq 'INDEL') {
#		$INFO{'Type'} = 'INDEL';
#		shift @groups;
#	} else {
#		$INFO{'Type'} = 'SNP';
#	}
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
	my (%GTitemCnt,$Mut,@plinkGT);
	for (@Samples) {
		if ($GT{$_}{'DP'} > 0 and $GT{$_}{'GQ'} >= 20) {
			my $gt = $GT{$_}{'GT'};
			++$GTcnt{$gt};
			my @GT = split /[\/|]/,$gt;
			++$SPcnt;
			if ($Pheno{$_} == 2) {
				++$GTitemCnt{$_} for @GT;
			}
			$GT{$_}{'GTp'} = join ' ',map($Bases[$_], @GT);
			#$GT{$_}{'O_K'} = 1;
		} else {
			$GT{$_}{'GTp'} = '0 0';
			#$GT{$_}{'O_K'} = 0;
		}
		push @plinkGT,$GT{$_}{'GTp'};
	}
	($Mut) = sort { $GTitemCnt{$b} <=> $GTitemCnt{$a} } keys %GTitemCnt;
	unless ($Mut) {
		++$Stat{'VCF_noMUT_Skip'};
		next;
	}
	$Mut = $Bases[$Mut];
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
	if ($QUAL<20 or $INFO{'FQ'}<=0 or scalar(keys %GTcnt)<2 or $SPcnt<6) {
		++$Stat{'VCF_Skipped'};
		next;
	}
	my $SNPid = "r".$Stat{'VCF_In'};
	$CHROM =~ /(\d+)/;
	print OP join("\t",$1,$SNPid,0,$POS,@plinkGT),"\n";
	print OM join("\t",$SNPid,$Mut),"\n";
	print OD join("\t",${CHROM},${POS},$SNPid),"\n";
}
close $th;

close OP;
close OM;
close OD;

ddx \%Stat;

print "[Prepare $outfs.phe.] And then:\np-link --tfile $outfs --reference-allele $outfs.MinorAllele --fisher --out ${outfs}P --model --cell 0\n--pheno $outfs.phe --all-pheno [screen log will be saved by p-link itselt]\n";
__END__
grep -hv \# radseq.gt > radseq.tfam

./bcf2bed.pl radseq.tfam radseq.bcgv.bcf radseq 2>&1 | tee radseq.pedlog

grep REC radseq.p.snow.model > radseq.p.snow.model.REC
grep DOM radseq.p.snow.model > radseq.p.snow.model.DOM
sort -nk8 radseq.p.snow.model.REC > radseq.p.snow.model.REC.sortnk8 &
sort -nk8 radseq.p.snow.model.DOM > radseq.p.snow.model.DOM.sortnk8 &
bcftools view -I radseq.bcgv.bcf |grep -v \# |cat -n > radseq.bcgv.bcf.rs &

bcftools view tigers.bcgv.bcf|grep -v \## > tigers.bcgv.vcf &
