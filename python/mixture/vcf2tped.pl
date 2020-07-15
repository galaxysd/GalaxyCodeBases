#!/usr/bin/env perl
use strict;
use warnings;

use Data::Dump qw(ddx);

my $Usage = "Usage: $0 <vcf.gz> <out prefix>\n";
die $Usage if @ARGV < 1;

my ($filename,$outp) = @ARGV;
my $cmd = "bcftools view -U -m2 -i '%QUAL>=40 & MIN(FMT/GQ)>20 & GT[2]=\"het\"' -v snps $filename |bcftools view -e 'FMT/DP=\".\"'| bcftools query -f '%CHROM\t%POS\t%REF,%ALT\t%QUAL[\t%TGT:%AD:%GQ]\n'";
my @SampleIDs;
open(SID,'-|',"bcftools query -l $filename") or die "Error opening [$filename]: $!\n";
while(<SID>) {
	chomp;
	push @SampleIDs,$_;
}
close SID;
ddx \@SampleIDs;

my $theKiller = 'T3C';
my $theVictim = 'T10C';
my @addedSID = qw(Rsin Rdip);
my @AllSID = (@SampleIDs,@addedSID);

open O,'>',"$outp.tfam" or die "[$outp.tfam]:$!\n";
for (@AllSID) {
	my $FID = 'T';
	my $IID = $_;
	print O join(' ',$FID,$IID,0,0,0,0),"\n";
	# https://www.cog-genomics.org/plink/1.9/formats#fam
}
close O;

my %COUNTING;
sub doCnt($) {
	my $s = $_[0];
	++$COUNTING{$s};
	print " $s";
}
sub mergeGT($) {
	my $hashref = $_[0];
	my @GTs = sort keys %{$hashref};
	if (scalar @GTs ==1) {
		push @GTs,@GTs;
	}
	return join(' ',@GTs);
}
open O,'>',"$outp.tped" or die "[$outp.tped]:$!\n";
open(IN,"-|",$cmd) or die "Error opening [$filename]: $!\n";
while (<IN>) {
	chomp;
	my ($Chrom,$Pos,$RefAlt,$Qual,@GTs) = split /\t/,$_;
	next if $Chrom =~ /_/;
	print "$_\n";
	my (%GTs,%GTstr,%GTped);
	@GTstr{@SampleIDs} = @GTs;
	#ddx \%GTstr;
	for my $k (keys %GTstr) {
		my ($sGT,$sAD,$sGQ) = split /:/,$GTstr{$k};
		my @aGT = split /[|\/]/,$sGT;
		my @aAD = split /\,/,$sAD;
		my %counter;
		for (@aGT) {
			++$counter{$_};
		}
		#my @result = sort(keys %counter);
		$GTs{$k} = \%counter;
	}
	$GTs{$addedSID[0]} = {%{$GTs{'mixed'}}};
	$GTs{$addedSID[1]} = {%{$GTs{'mixed'}}};
	for (keys %{$GTs{$theVictim}}) {
		if (exists $GTs{$addedSID[0]}{$_}) {
			delete $GTs{$addedSID[0]}{$_};
		}
	}
	for my $k (keys %GTs) {
		$GTped{$k} = mergeGT($GTs{$k});
	}
	if (length($GTped{$addedSID[0]})==0) {
		doCnt('0same');
		next;
	}
	ddx \%GTs;
	ddx \%GTped;
}
close IN;
# https://www.cog-genomics.org/plink/1.9/formats#tped
close O;
