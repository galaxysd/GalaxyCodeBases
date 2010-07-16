#!/bin/env perl
use strict;
use warnings;
#use Data::Dump qw(dump ddx);

unless (@ARGV){
	print "perl $0 <chr.nfo> <soaps.stat> <out file>\n";	# soaps.nfo can store the file size of soap. Too late, useless.
	exit;
}

my ($nfo,$stat,$out) = @ARGV;

my (%ChrLen,%Dat,@Samples);
open N,'<',$nfo or die "[!]Error opening $nfo: $!\n";
while(<N>) {
	my ($chr,$len)=split /\t/;
	$ChrLen{$chr}=0+$len;
}
close N;
#ddx \%ChrLen;

open N,'<',$stat or die "[!]Error opening $stat: $!\n";
my %samples;
while(<N>) {
	next if /^(ALL|#)/;
	next if /^\s*$/;
	my ($chr,$sample,$dat)=split /\t/;
	$Dat{$chr}{$sample}=(split /,/,$dat)[1] if $dat=~/\d/;
	++$samples{$sample};
}
close N;
@Samples = sort keys %samples;
#ddx \%Dat;

open O,'>',$out or die "[!]Error opening $out: $!\n";
print O join("\t",'ChrID','ChrLength',@Samples),"\n";
for my $chr (sort { $ChrLen{$b} <=> $ChrLen{$a} } keys %Dat) {
	my $len=$ChrLen{$chr};
	print O $chr,"\t",$len;
	for my $sample (@Samples) {
		my $t=$Dat{$chr}{$sample} || 0;
		print O "\t",$t,',',sprintf("%.2f",$t/$len);
	}
	print O "\n";
}
close O;
