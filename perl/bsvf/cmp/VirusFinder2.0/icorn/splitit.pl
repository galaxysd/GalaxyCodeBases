#! /usr/bin/perl -w
#
# File: splitit.pl
# Time-stamp: <22-Oct-2009 17:13:40 tdo>
# $Id: $
#
# Copyright (C) 2008 by Pathogene Group, Sanger Center
#
# Author: Thomas Dan Otto
#
# Description: Will spilt a given file in the pieces with fixed amount
# of lines and read them with prefix resultName.
#

use strict;

if (scalar(@ARGV)<3) {
  die "Usage: <fastq/a file> <amount of lines> <resultName>\n"
}

$re= Helps::splitit(shift,shift,shift);

print $re;

sub splitit{
  
  my $name       = shift;
  my $amount     = shift;
  my $resultname = shift;
  
  open (F, $name) or die "couldn't open file $name: $! \n";

  my $count=0;
  my $res='';
  my $fileCount=1;
  ## need clean fastq format, no return
  while (<F>) {
	$res.=$_;
	$res.=<F>;
	$res.=<F>;
	$res.=<F>;
	$count++;
	if (($count%$amount)==0) {
	  open (WI, "> $resultname.$fileCount") or die "couldn't write...\n";
	  print WI $res;
	  close(WI);
	  $fileCount++;
	  
	  $res=''; 
	}
  }
  
  if ($res ne '') {
	open (WI, "> $resultname.$fileCount") or die "couldn't write...\n";
	print WI $res;
	close(WI);
	
	return $fileCount;
  }
  else {
	return ($fileCount-1);
  }
}



