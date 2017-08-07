#! /usr/bin/perl -w
#
# File: join.solexa2pairs.pl
# Time-stamp: <03-Jul-2009 12:19:53 tdo>
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


my $name1=shift;
my $resultName = shift;

if (!defined($resultName)) {
  die "Parameter: <root_s_x.mate.fastq> <Resultname>\n";
  
}


open (F1,$name1) or die "Problems to open $name1";

my $count=0;
my $res;

open (F, "> $resultName" ) or die "probs";
my $seqtmp;
my $qualtmp;

while (<F1>) {
  ### change /1 to .F
  s/\.F$//g;

  # name is the same
  $res .= $_;
  $count++;

  
  $seqtmp =<F1>;
  chomp($seqtmp);
  $_=<F1>;
  $qualtmp =<F1>;
  chomp($qualtmp);

  #line of reverse
  $_=<F1>;

  #seq
  $res .=$seqtmp.<F1>;
  $count++;
  $res.=<F1>;

  $res .=$qualtmp.<F1>;
  
  if (($count%100000)==0) {
	print F $res;
	$res='';
	
  }
}

print F $res;
close(F);

