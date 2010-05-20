#!/bin/env perl
use strict;
use warnings;

unless (@ARGV){
	print "perl $0 <fq pe list> <fq se list> <fq stats>\n";
	exit;
}

my ($pelst,$selst,$statout) = @ARGV;
my (%DATrbrf,%Librmx,%FQnfo);
open LST,'<',$pelst or die "Error opening $pelst: $!\n";
while (<LST>) {
	chomp;
	my ($sample,$lib,$FL,$min,$max,$fq1,$fq2,$ext,$path)=split /\t/;
	$Librmx{$sample}{$lib}{$FL}=[0,$min,$max];	# readlen,minIS,maxIS
	$FQnfo{$sample}{$lib}{$FL}=[$ext,$path,$fq1,$fq2];
	$DATrbrf{$sample}{$lib}=[0,0,0,0];	# rawReads,rawBP,filteredReads,filteredBP
}
open LST,'<',$selst or warn "Error opening $pelst: $!\n";
while (<LST>) {
	chomp;
	my ($sample,$lib,$FL,$fq,$ext,$path)=split /\t/;
	$Librmx{$sample}{$lib}{$FL}=[0];	# readlen,minIS,maxIS
	$FQnfo{$sample}{$lib}{$FL}=[$ext,$path,$fq];
	$DATrbrf{$sample}{$lib}=[0,0,0,0];	# rawReads,rawBP,filteredReads,filteredBP
}

for my $sample (keys %FQnfo) {
	for my $lib (keys %{$FQnfo{$sample}}) {
		for my $FL (keys %{$FQnfo{$sample}{$lib}}) {
			my ($ext,$path,@fq)=@{$FQnfo{$sample}{$lib}{$FL}};
			for (@fq) {
				;
			}
		}
	}
}