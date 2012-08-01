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
use List::MoreUtils 'first_index';

die "Usage: $0 <plink model file> <bcgv bcf> <out>\n" if @ARGV<3;
my $plinkfs=shift;
my $bcfs=shift;
my $outfs=shift;

my (%Stat,$t);
open O,'>',$outfs or die "Error opening $outfs: $!\n";
$t = "# In: [$plinkfs],[$bcfs], Out: [$outfs]\n";
print O $t;
print $t;

my (%Plink,@PlinkS,%Vcf,$SortI);
open L,'<',$plinkfs or die;
while (<L>) {
	if (/^ CHR /) {
		s/^\s+//;
		$SortI = first_index { $_ eq 'P' } (split /\s+/);
#warn $SortI,(split /\s+/)[$SortI];
		next;
	}
	chomp;
	s/^\s+//;
	my @Dat = split /\s+/;	# $chr,$markID,$min,$maj,$model,@other
	$Dat[1] =~ /(\d+)/;
	$Dat[1] = $1;
#die '[',join('|',@Dat),"]\n";
	$Plink{$Dat[1]} = \@Dat;
}
close L;
@PlinkS = sort { $Plink{$a}->[$t] <=> $Plink{$b}->[$t] } keys %Plink;

my $th = openpipe('bcftools view -I',$bcfs);	# -I	skip indels
my (@Samples,@Parents);
while (<$th>) {
	next if /^##/;
	chomp;
	my @data = split /\t/;
	if ($data[0] eq '#CHROM') {
		@Samples = map {my $t=(split /\//)[-1];$t=~s/_cut//g;$t=~s/-/./g; $_=join('_',(split /\./,$t)[-5,-6]);} splice @data,9;
		# ../5.bam_0000210210_merged/d1_4_merged.D4.JHH001.XTU.sort.rmdup.bam
		@Parents = grep(!/^GZXJ/,@Samples);
		last;
	}
}
print O "# Samples: [",join('],[',@Samples),"]\n# Parents: [",join('],[',@Parents),"]\n";
warn "Samples:\n[",join("]\n[",@Samples),"]\nParents: [",join('],[',@Parents),"]\n";

my $VCF_In;
while (<$th>) {
	next if /^#/;
	chomp;
	#my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, @data) = split /\t/;
	++$VCF_In;	# also as rs#
	if (exists $Plink{$VCF_In}) {
		$Vcf{$VCF_In} = $_;
	}
}
close $th;
close O;
