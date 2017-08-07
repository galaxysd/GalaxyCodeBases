#! /usr/bin/perl -w
#
# File: analyse.roundMatch.pl
# Time-stamp: <15-May-2009 13:03:01 tdo>
# $Id: $
#
# Copyright (C) 2008 by Pathogene Group, Sanger Center
#
# Author: Thomas Dan Otto
#
# Description:
#


my $cigar = shift;


open (F, $cigar) or die "proble\n";

my @ar=<F>;

my $LIM=shift;
my $minInsertSize=shift;

if (!defined($LIM)) {
  $LIM=1000;
}
if (!defined($minInsertSize)) {
  $minInsertSize=50;
}
my %h;
my $h2;
my %Dist;
my $goodPair;

my $count=0;

my $dist;
my $basesMatch=0;


foreach (@ar) {
  my @ar=split(/\s+/);
  my ($root) = $ar[1] =~ /^(.*)\./;
  $h{$root}++;
  $basesMatch+=$ar[9];
  
  if (!defined($Dist{$root})) {
	$Dist{$root}{"chr"}=$ar[5];
	$Dist{$root}{"strand"}=$ar[4];
	if ($ar[4] eq '-') {
	  $Dist{$root}{"pos"}=$ar[7];
	}
	else
	  {
		$Dist{$root}{"pos"}=$ar[6];
	  }
	
  }
  else {
	$count++;
	my $maxHitRef=$ar[6];
	if ($ar[4] eq '-') {
	  $maxHitRef=$ar[7];
	}
	my $insertSize=(abs($Dist{$root}{"pos"} - $maxHitRef));
	
	if ( ($Dist{$root}{"chr"} eq $ar[5]) and
		 ($Dist{$root}{"strand"} ne $ar[4]) and
		 ( $insertSize > $minInsertSize ) and
		 ( $insertSize < $LIM )
	   )
	  {
		$goodPair++;
		$dist+= ($insertSize);
		
	  }
  }
}




print "$count\n$goodPair\n";

printf ("%.0f\n", ($dist/$goodPair));
print "$basesMatch\n";


sub getAnchor{
  my $name=shift;


  open (F, $name) or die "couldn open anchor cigar file $!\n";

  my @ar=<F>;

  my %h;
  
  foreach (@ar) {
	my @a =split(/\s/);
	$h{$a[1]}{pos}=$a[6];
	$h{$a[1]}{chr}=$a[5];
	$h{$a[1]}{strand}=$a[4];
	
  }  
  return(\%h);
}
