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

my (%Plink,@PlinkS,%Vcf,$SortI,$TypeI,%Types);
open L,'<',$plinkfs or die;
while (<L>) {
	if (/^ CHR /) {
		s/^\s+//;
		$SortI = first_index { $_ eq 'P' } (split /\s+/);
		$TypeI = first_index { $_ eq 'TEST' } (split /\s+/);
#warn $SortI,(split /\s+/)[$SortI];
		next;
	}
	chomp;
	s/^\s+//;
	my @Dat = split /\s+/;	# $chr,$markID,$min,$maj,$model,@other
	$Dat[1] =~ /(\d+)/;
	$t = $1;
#die '[',join('|',@Dat),"]\n";
	$Plink{$t}{$Dat[$TypeI]} = \@Dat;
	++$Types{$TypeI};
}
close L;
ddx \%Plink;
die;

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

my ($VCF_In,%ChrPn,%ChrPs,%TypePsum);
for (keys %Types) {
	$ChrPs{$_} = ();
}
while (<$th>) {
	next if /^#/;
	chomp;
	#my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, @data) = split /\t/;
	my @data = split /\t/;
	++$VCF_In;	# also as rs#
	if (exists $Plink{$VCF_In}) {
		$Vcf{$VCF_In} = \@data;
		for (keys %{$Plink{$VCF_In}{$_}}) {
			$t = $Plink{$VCF_In}{$_}->[$SortI];
			$ChrPs{$_}{$data[0]} += $t;
			$TypePsum{$_} += $t;
		}
		++$ChrPn{$data[0]};
		#$ChrPs{$data[0]} += $Plink{$VCF_In}->[$SortI];
	}
}
close $th;
warn "bcf done.\n";

print O '# ChrID Count: ',scalar(keys %ChrPn),"\n",'# SNP Count: ',scalar(keys %Plink),"\n";

for my $type (sort { $TypePsum{$a} <=> $TypePsum{$b} } keys %Types) {
	for my $chr (keys %ChrPn) {
		$ChrPs{$type}{$chr} /= $ChrPn{$chr};
	}
	@PlinkS = sort { $ChrPs{$type}{$Vcf{$a}->[0]} cmp $ChrPs{$type}{$Vcf{$b}->[0]} || $Plink{$a}{$type}->[$SortI] <=> $Plink{$b}{$type}->[$SortI] } keys %Plink;
	for (@PlinkS) {
		print O join("\t",@{$Plink{$_}{$type}},@{$Vcf{$_}}),"\n";
	}
}
close O;
