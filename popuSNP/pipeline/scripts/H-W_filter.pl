#!/bin/env perl
use strict;
use warnings;
use Getopt::Long;

#######################
##USAGE
#######################

sub usage{
	print STDERR <<USAGE;
Version:1.0
2009-10-12       hewm\@genomics.org.cn
ReVersion:1.01
Mon Oct 12 14:09:56 CST 2009

		Usage: perl $0
		where:

		Options
			-snp  <c>: (must be gave)
			-rm   <c>: (must be gave)
			-o    <c>: outDir
			-h       : show this help message
USAGE
}

my %opt;
GetOptions(\%opt,"snp=s","o=s","rm=s","h");

#########check parameters and write into memory ######
if(exists $opt{h}){
	usage;
	exit;
}
unless($opt{snp} || $opt{rm}){
	usage;
	print STDERR "\n\t\t\tERROR: -i option missing!\n";
	exit;
}

my $temp=`pwd`;  chomp $temp ; $temp=$temp."\/out.out";


$opt{o} ||=$temp;

##############what you want to do #################
open A , "$opt{rm}" || die "the Input Can't open $!";

my %hash=();

 while(<A>)
 {
     chomp ;
     my @inf=split ;
     $hash{$inf[1]}=1;
 }
 close A;

open A ,'<', $opt{snp} || die "the Input Can't open $!";
open B ,'>', $opt{o} || die "the Input Can't open $!";

 while(<A>)
 {   my $line=$_ ;
     chomp ;
     my @inf=split ;
     if(!exists $hash{$inf[1]})
     {
         print B  $line ;
     }
 }

 close B ;
 close A ;
### swimming in the sky and flying in the sea ####
