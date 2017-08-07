#! /usr/bin/perl -w
#
# File: join.solexa2pairs.pl
# Time-stamp: <03-Jul-2009 12:19:16 tdo>
# $Id: $
#
# Copyright (C) 2008 by Pathogene Group, Sanger Center
#
# Author: Thomas Dan Otto
#
# Description: Will get the two solexa files \1 and \2 and join them,
# replacing the name to .F and .R of each read name
#

use strict;
use warnings;

my $name1=shift;
my $name2=shift;
my $resultName = shift;

print "Transform the two fastq file to the correct fileformat for SSAHA-pileup\n";
if (!defined($resultName)) {
  die "Parameter: <forward solexa file (\1)> <reverse solexa file (\2)> <Resultname>\n";
  
}

if ($name2 =~ /\.gz/) {
  open F2," gunzip -c $name2 | " or die "cannot open $name2\n";
}
else {
  open F2,$name2 or die "cannot open $name2\n";
}
my @ar=<F2>;
close(F2);

if ($name1 =~ /\.gz/) {
        open F1," gunzip -c $name1 | " or die "cannot open $name1\n";
  }
  else {
        open F1,$name1 or die "cannot open $name1\n";
  }


my $count=0;
my $res;

open (F, "> $resultName" ) or die "problems to write the output file $resultName: $!\n";

while (<F1>) {
  ### change /1 to .F
  if (/\/1$/){
	  s/\/1$/\.F/g;
	  $res .= $_;
  }
  else
  { 
  	chomp;
	$res .= $_.".F\n";
  }
  
  $res .=<F1>; $res .=<F1>; $res .=<F1>;

  ## get the reverse
  if ($ar[$count] =~ /\/2$/){
	 $ar[$count] =~ s/\/2$/\.R/g;
  }
  else 
  {
  	chomp($ar[$count]);
	$ar[$count] .= ".R\n";
  }
  
  $res .=$ar[$count++];
  $res .=$ar[$count++];$res .=$ar[$count++];$res .=$ar[$count++];

  if (($count%10000)==0 ||($count%10000)==1 ||($count%10000)==2 ||($count%10000)==3 ) {
	print F $res;
	$res='';
	
  }
}

print F $res;
close(F);

