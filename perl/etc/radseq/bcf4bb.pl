#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Data::Dump qw(dump ddx);
use Galaxy::IO;
use Galaxy::SeqTools;

die "Usage: $0 <bcgv bcf> <out>\n" if @ARGV<2;
my $bcfs=shift;
my $outfs=shift;

my (%Stat,$t);
open O,'>',$outfs.'.out' or die "Error opening $outfs.out : $!\n";
$t = "# In: [$bcfs], Out: [$outfs]\n";
print O $t;
print $t;

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
			++$GTitemCnt{$_} for @GT;
			$GT{$_}{'GTp'} = join ' ',map($Bases[$_], @GT);
			#$GT{$_}{'O_K'} = 1;
		} else {
			$GT{$_}{'GTp'} = '0 0';
			#$GT{$_}{'O_K'} = 0;
		}
		push @plinkGT,$GT{$_}{'GTp'};
	}
	($Mut) = sort { $GTitemCnt{$b} <=> $GTitemCnt{$a} } keys %GTitemCnt;
	next if (keys %GTitemCnt) > 1;
	unless (defined $Mut) {
		++$Stat{'VCF_noAffInd_Skipped'};
		next;
	}
	unless ($Mut) {
		++$Stat{'VCF_noMUT_Skipped'};
		next;
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
	if ($QUAL<20 or $INFO{'FQ'}<=0 or $SPcnt<17) {
		++$Stat{'VCF_Skipped'};
		next;
	}
	$Mut = $Bases[$Mut];
	my $SNPid = "r".$Stat{'VCF_In'};
	$CHROM =~ /(\d+)/;
	print STDERR "[$SPcnt] $CHROM, $POS, $REF, $ALT, $QUAL, $INFO ",dump(\%GT),"\n";
}
close $th;

ddx \%Stat;

close O;
__END__

./bcf4bb.pl tigers.bcgv.bcf t 2>&1|tee t.log

grep scaffold b17.log |awk '{print $1}'|uniq|perl -lane 's/\,//;print $_'|while read a;do cat bychr/$a.fa >> bgene.fa;done

cat bgene.chr|while read a;do cat P_tigris.gene.gtf |grep -P "$a\t" >> bgene.vcf ;done
