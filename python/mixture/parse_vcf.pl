#!/usr/bin/env perl
use strict;
use warnings;

use Data::Dump qw(ddx);

my $Usage = "Usage: $0 <vcf.gz>\n";
die $Usage if @ARGV < 1;

my ($filename) = @ARGV;
#my $cmd = "bcftools view -m3 -v snps $filename | bcftools query -f '%CHROM\t%POS\t%REF,%ALT\t%QUAL[\t%TGT:%AD:%GQ]\n'";
my $cmd = "bcftools view -U -m2 -i '%QUAL>=40 & MIN(FMT/GQ)>20 & GT[2]=\"het\"' -v snps $filename |bcftools view -e 'FMT/DP=\".\"'| bcftools query -f '%CHROM\t%POS\t%REF,%ALT\t%QUAL[\t%TGT:%AD:%GQ]\n'";
=pod
chr1	102951	C,T	1575.58	C/C/C/T:23,6:25	C/C/C/T:88,27:67	C/C/C/T:111,33:93
chr1	633963	C,T	6936.96	T/T/T/T:5,242:99	C/C/C/C:250,0:99	C/C/C/T:172,77:75
chr1	889636	A,G,T	1356.65	A/A/G/T:14,3,8:21	A/A/G/T:36,20,10:48	A/A/G/T:50,23,18:75
chr1	1099396	G,A,C	214.66	G/G/G/G:22,1,0:24	G/G/A/C:67,5,9:58	G/G/A/C:89,6,9:24
=cut
my @SampleIDs;
open(SID,'-|',"bcftools query -l $filename") or die "Error opening [$filename]: $!\n";
while(<SID>) {
	chomp;
	push @SampleIDs,$_;
}
close SID;
ddx \@SampleIDs;
# ["T10C", "T3C", "mixed"] => v,k,m

sub mergeGT($) {
	my $hashref = $_[0];
	my @GTs = sort keys %{$hashref};
	return join('',@GTs);
}

my %COUNTING;
open(IN,"-|",$cmd) or die "Error opening [$filename]: $!\n";
while (<IN>) {
	chomp;
	my ($Chrom,$Pos,$RefAlt,$Qual,@GTs) = split /\t/,$_;
	next if $Chrom =~ /_/;
	print "$_\n";
	my %GTs;
	for my $i (0 .. $#SampleIDs) {
		my ($sGT,$sAD,$sGQ) = split /:/,$GTs[$i];
		my @aGT = split /[|\/]/,$sGT;
		my @aAD = split /\,/,$sAD;
		my %counter;
		for (@aGT) {
			++$counter{$_};
		}
		#my @result = sort(keys %counter);
		$GTs{$SampleIDs[$i]} = \%counter;
	}
	$GTs{'__R'} = {%{$GTs{'mixed'}}};
	for (keys %{$GTs{'T10C'}}) {
		if (exists $GTs{'__R'}{$_}) {
			delete $GTs{'__R'}{$_};
		}
	}
	$GTs{'_V'} = mergeGT($GTs{'T10C'});
	$GTs{'_K'} = mergeGT($GTs{'T3C'});
	$GTs{'_M'} = mergeGT($GTs{'mixed'});
	$GTs{'_R'} = mergeGT($GTs{'__R'});
	ddx \%GTs;
	++$COUNTING{'_All'};
	if ($GTs{'_R'} eq '') {
		++$COUNTING{'0noR'};
	} elsif ($GTs{'_R'} eq $GTs{'_K'}) {
		++$COUNTING{'1fulEQ'};
	} elsif (length($GTs{'_R'})==1 and $GTs{'_K'} =~ /$GTs{'_R'}/) {
		++$COUNTING{'2inc'};
	} else {
		++$COUNTING{'3ne'};
	}
}
ddx \%COUNTING;
close IN;
