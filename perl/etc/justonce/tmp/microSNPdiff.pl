#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

my $SampleList = 'samples.lst';
my $SNPtsv = 'micro.snp.tsv';

my (%Samples,@Samples,%ID2Sample,%SampleCnts);
open I,'<',$SampleList or die $?;
while (<I>) {
	chomp;
	my ($id,$rep) = /([A-Za-z]+)(\w+)/;
	push @Samples,$_;
	push @{$Samples{$1}},$_;
	$ID2Sample{$_} = $1;
	++$SampleCnts{$1};
}
close I;

my ();
open I,'<',$SNPtsv or die $?;
while (<I>) {
	chomp;
	my ($Chr,$Pos,$GTstr,$Qual,@Dat) = split /\t/;
	my @GTs = split /,/,$GTstr;
	my @DatA = map { (split /;/,$_)[0] } @Dat;
	#ddx [$Chr,$Pos,\@GTs,$Qual,\@DatA];
	my %SmpGTp;
	for my $i (0 .. $#Samples) {
		my $sid = $Samples[$i];
		my $uid = $ID2Sample{$sid};
		my $GT = $DatA[$i];
		++$SmpGTp{$uid}{$GT};
	}
	ddx \%SmpGTp;
}
close I;



__END__

ddx \%Samples;
ddx \%ID2Sample;
ddx \%SampleCnts;
print scalar keys %SampleCnts,"\n";

