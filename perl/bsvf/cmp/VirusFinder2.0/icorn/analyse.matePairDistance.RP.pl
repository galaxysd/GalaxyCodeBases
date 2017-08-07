#! /usr/bin/perl -w
#
# File: analyse.roundMatch.pl
# Time-stamp: <26-May-2009 16:02:54 tdo>
# $Id: $
#
# Copyright (C) 2008 by Pathogene Group, Sanger Center
#
# Author: Thomas Dan Otto
#
# Description:
#
# use the most filtered optino please!

use strict;

# maximal distance 
 
if (scalar(@ARGV)<3) {
  die "perl analyse.matePairDistance.pl <cigar> <Resultname> <MAX>\n";
}
my $cigar = shift;
my $resultName = shift;
my $MAX=shift;
my $MIN = shift;

if (!defined($MIN)){
   $MIN=50;
}

print "vai\n";

open (F, $cigar) or die "proble\n";

#my @F=<F>;

my %h;
my %h2;

my %res;
my $dist=0;


my $oldRoot='';

#we can assume, that the order is alrway root.p root.q nrew
my $count=0;
my $count1=0;
my $result;


my $lastline=<F>;


#for (my $line=0;$line<scalar(@F);$line++) {
while  (<F>) {
  my $line =$_;#  print $line ."\n";
  

  my ($root)  = $lastline   =~ /^cigar::\d+\s(\S+)\./;
  my ($rootR) = $line  =~ /^cigar::\d+\s(\S+)\./;
  #$arR[1] =~ /^(.*)\./;
  
  if ($root eq $rootR) {
	checkReadPair($lastline,$_);
	#	checkReadPair($F[($line)],$F[($line+1)]);
	#	$line++;
#	$_=<F>;
  }
  $lastline=$line;
  
}


if (0) {

my $res;

foreach my $chr (sort keys %res) {
  foreach my $pos (sort {$a <=> $b} keys %{ $res{$chr} }) {
	foreach my $read (keys %{ $res{$chr}{$pos} }) {
	  $res .= "$chr\t$pos\t$read\t$res{$chr}{$pos}{$read}";	  
	}	
  }
}
}


print "$count pairs were found\n$count1 pairs were between $MIN and $MAX bp \nMean is ";
printf ("%.0f ", ($dist/$count1));
print ("(just considering pairs dist > $MIN and < $MAX)\n");
#print "if awk '{ if (\$5> $MIN && \$5 < $MAX){ print \$5 }}' $resultName.allDistance.txt > forR.txt you get a list to do a plot with R (dat<-read.table(\"forR.txt\");plot(density(dat[,])); , or so \n";

open (F, "> $resultName.cigar") or die "please give resultnam $resultName \n";
print F $result;
close(F);



#open (F, "> $resultName.allDistance.txt") or die "please give resultnam $resultName \n";
#print F $res;
#close(F);



sub checkReadPair{
  
  my $p=shift;
  my $q=shift;

  my ($chr1,$pos1s,$pos1e,$strand1,$matches1,$read1start,$read1end)=getPos( $p);
  my ($chr2,$pos2s,$pos2e,$strand2,$matches2,$read2start,$read2end)=getPos( $q);
  
  $count++;
  if ($strand1 ne $strand2 
	  and
	  $chr1 eq $chr2
	  and
	  (
	   ($strand1 eq '+' and $pos1s < $pos2s # read is left
		and $read1start < $read1end   # not complement positive
		and $read2start > $read2end
	   )
	   or
	   ($strand2 eq '+' and $pos2s < $pos1s # read is left
		and $read2start < $read2end   # not complement positive
		and $read1start > $read1end 
	   )
	  )
	 ) {
	my @sorted = sort { $a <=> $b} ($pos1s,$pos1e,$pos2s,$pos2e);
	#	$res{$chr1}{$sorted[0]}{$p}.="$q\t".($sorted[3]-$sorted[0])."\t$matches1\t$matches2\n";
	
	if (abs($sorted[3]-$sorted[0])>$MIN
		and
		abs($sorted[3]-$sorted[0])<$MAX
	   ) {
	  
	  $dist+=($sorted[3]-$sorted[0]);
	  $count1++;
	  $result.=$p;
	  $result.=$q;
	  #	  print "$p $q";
	  
	}
	
  }
  else {
#	print "$p$q";
	
  }
}

sub getPos{
  my $line = shift;

  chomp($line);
  my @ar=split(/\s/,$line);

  return ($ar[5],$ar[6],$ar[7],$ar[4],$ar[9],$ar[2],$ar[3]);
  


}
