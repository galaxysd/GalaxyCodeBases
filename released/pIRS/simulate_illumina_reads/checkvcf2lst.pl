#!/bin/env perl
use strict;
use warnings;

unless (@ARGV){
	print "Usage $0 <sim-list-prefix> <vcf>\n";
	exit;
}

my (%SNP,%Indel);
my ($ccSNP,$ccIndel)=(0,0);

open IS,'<',$ARGV[0].'_snp.lst' or die "Error: $!\n";
while (<IS>) {
	chomp;
	my ($chr,$pos,$ref,$snp)=split /\t/;
	$SNP{$chr}{$pos}=$ref."\t".$snp;
	++$ccSNP;
}
close IS;
open INDEL,'<',$ARGV[0].'_indel.lst' or die "Error: $!\n";
while (<INDEL>) {
	chomp;
	my ($chr,$pos,$insdel,$count,$seq)=split /\t/;
	$Indel{$chr}{$pos}=join("\t",$insdel,$count,$seq);
	++$ccIndel;
}
close INDEL;

my @cSNP=(0,0,0,0,0);	# Error TP TN FP FN
my @cIndel=(0,0,0,0,0);

sub vcf2indel($$) {
	my ($ref,$alt)=@_;
	my $lref=length $ref;
	my $lalt=length $alt;
	my ($insdel,$count,$seq);
	if ($lref>$lalt) {
		$insdel='-';
		$count = $lref - $lalt;
		$seq=substr $ref,$lalt,$count;
	} else {
		$insdel='+';
		$count = $lalt - $lref;
		$seq=substr $alt,$lref,$count;
	}
	return join("\t",$insdel,$count,$seq);
}

open VCF,'<',$ARGV[1] or die "Error: $!\n";
open DUP,'>',$ARGV[1].'.dup' or die "Error: $!\n";
while (<VCF>) {
	next if /^#/;
	my ($chr,$pos,$ref,$alt,$Q,$info)=(split /\t/)[0,1,3,4,5,7];
	if ($info =~ /\bINDEL\b/) {
		#print;\
		unless (exists $Indel{$chr}{$pos}) {
			++$cIndel[3];
			print DUP "FPIndel: $_";
			next;
		}
		#my ($insdel,$count,$seq)=split /\t/,$Indel{$chr}{$pos};
		my $ret=&vcf2indel($ref,$alt);
		if ($ret eq $Indel{$chr}{$pos}) { ++$cIndel[1]; }
		else { ++$cIndel[0]; print DUP "ErrIndel: [$Indel{$chr}{$pos}] <$ret>:$_"; }
	} else {
		unless (exists $SNP{$chr}{$pos}) {
			++$cSNP[3];
			print DUP "FPSNP: $_";
			next;
		}
		#my ($Sref,$Ssnp)=split /\t/,$SNP{$chr}{$pos};
		if ( $ref."\t".$alt eq $SNP{$chr}{$pos} ) { ++$cSNP[1]; }
		else { ++$cSNP[0]; print DUP "ErrSNP: [$SNP{$chr}{$pos}] $_"; }
	}
}
close VCF;
$cSNP[4] = $ccSNP - $cSNP[1] - $cSNP[0];
$cIndel[4] = $ccIndel - $cIndel[1] - $cIndel[0];
print "SNP   Error,TP,TN,FP,FN = ",join(", ",@cSNP),"\n";
print "Indel Error,TP,TN,FP,FN = ",join(", ",@cIndel),"\n";
print "---Indel results may be wrong when alt-seq forms direct repeat with terminals---\n";
close DUP;

