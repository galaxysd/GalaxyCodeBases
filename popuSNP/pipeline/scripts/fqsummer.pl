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
	$DATrbrf{Lib}{$sample}{$lib}=[0,0,0,0];	# rawReads,rawBP,filteredReads,filteredBP
	$DATrbrf{Sample}{$sample}=[0,0,0,0];
	$DATrbrf{Lane}{$sample}{$lib}{$FL}=[0,0,0,0];
}
open LST,'<',$selst or warn "Error opening $pelst: $!\n";
while (<LST>) {
	chomp;
	my ($sample,$lib,$FL,$fq,$ext,$path)=split /\t/;
	$Librmx{$sample}{$lib}{$FL}=[0];	# readlen,minIS,maxIS
	$FQnfo{$sample}{$lib}{$FL}=[$ext,$path,$fq];
	$DATrbrf{Lib}{$sample}{$lib}=[0,0,0,0];	# rawReads,rawBP,filteredReads,filteredBP
	$DATrbrf{Sample}{$sample}=[0,0,0,0];
	$DATrbrf{Lane}{$sample}{$lib}{$FL}=[0,0,0,0];
}

sub sumup ($$) {
	my ($arrayr,$hashr)=@_;
	$$arrayr[0] += $$hashr{InReads};
	$$arrayr[1] += $$hashr{InBPs};
	$$arrayr[2] += $$hashr{OutReads};
	$$arrayr[3] += $$hashr{OutBP};
}

my (@NFO,%NFO,$arrayp);
open OP,'>',$pelst.'.n' or die "Error opening $pelst.n: $!\n";
open OS,'>',$selst.'.n' or die "Error opening $selst.n: $!\n";
for my $sample (keys %FQnfo) {
	for my $lib (keys %{$FQnfo{$sample}}) {
		for my $FL (keys %{$FQnfo{$sample}{$lib}}) {
			my ($ext,$path,@fq)=@{$FQnfo{$sample}{$lib}{$FL}};
			for (@fq) {
				#print "[$_]\n";
				open NFO,'<',"$path$_.nfo" or die "Error opening $path$_.nfo: $!\n";
				@NFO=<NFO>;
				chomp @NFO;
				%NFO = map {split /\t/} grep /\t/,@NFO;
				#print "$_ => $NFO{$_}\n" for keys %NFO;
				close NFO;
				$Librmx{$sample}{$lib}{$FL}->[0]=$NFO{MaxReadLen} if $Librmx{$sample}{$lib}{$FL}->[0] < $NFO{MaxReadLen};
				&sumup($DATrbrf{Sample}{$sample},\%NFO);
				&sumup($DATrbrf{Lib}{$sample}{$lib},\%NFO);
				&sumup($DATrbrf{Lane}{$sample}{$lib}{$FL},\%NFO);
			}
			if ($#fq == 1) {	# PE
				print OP join("\t",$sample,$lib,$FL,@{$Librmx{$sample}{$lib}{$FL}},@{$FQnfo{$sample}{$lib}{$FL}}),"\n";
			} else {
				print OS join("\t",$sample,$lib,$FL,@{$Librmx{$sample}{$lib}{$FL}},@{$FQnfo{$sample}{$lib}{$FL}}),"\n";
			}
		}
	}
}
close OP; close OS;

open O,'>',$statout or die "Error opening $statout: $!\n";

close O;
