#!/bin/env perl
use strict;
use warnings;
# original from HeWeiMing;   Wed Oct 14 15:16:22 CST 2009

die"Version 1.1\t2010-7-22;\nUsage: $0 <GenotyepSnp> <Ref_chr.fa> <Out file>\n" unless (@ARGV ==3);

#############Befor  Start  , open the files ####################

open IA,'<',$ARGV[0]  || die "input file can't open $!" ;
open IB,'<',$ARGV[1]  || die "input file can't open $!" ;
open OA,'>',$ARGV[2] || die "output file can't open $!" ;

################ Do what you want to do #######################
$/=">";
<IB> ;
my $ref=();
while(<IB>)
{
	chomp ;
	my @inf=split /\n/;
	$ref=join("",@inf[1..$#inf]) ;
}
$/="\n";

my $header=<IA>;
if ($header =~ /^#ChrID/) {
	my $Samples=(split /\t/,$header)[2];
	$header="#ChrID\tPos\tRef\t$Samples";	# Well, no chomp no "\n".
} else {
	seek IA,0,0;
	$header=undef;
}
if (defined $header) {
	print OA $header;
}

while(<IA>)
{
	chomp ;
	my @inf=split /\t/;
	my $lie=substr($ref,$inf[1]-1,1);
	print OA $inf[0],"\t",$inf[1],"\t",$lie,"\t",$inf[2],"\n";
}

close OA ;
close IB;
close IA;
##############swiming in the sky and flying in the sea #######################
