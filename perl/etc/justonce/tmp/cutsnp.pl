#!/usr/bin/env perl
use strict;
use warnings;
use Text::CSV;
use Data::Dump qw(ddx);

my $dbSNPa = 'snp151.lst.h';

my %rsIN;
open IN,'<','bgi2.tsv' or die "Error opening file: $!\n";
<IN>;
while(<IN>) {
	chomp;
	my ($rsid,$chr,$Ref,$Alt,$pGlobal,$pAsia) = split /\s+/;
	$rsIN{$rsid} = [$chr,$Ref,$Alt,$pGlobal,$pAsia];
}
close IN; 
ddx \%rsIN;

my %rsDB;
open R,'<',$dbSNPa or die "Error opening [$dbSNPa]: $!\n";
while(<R>) {
	chomp;
	my ($chr,$pos,$rsid,$ref,$tgt) = split /\s+/;
	if (exists $rsIN{$rsid}) {
		my @dat = @{$rsIN{$rsid}};
		$rsDB{$rsid} = [$chr,$pos,$ref,$tgt,@dat[-2,-1]];
	}
}
close R;
ddx \%rsDB;
