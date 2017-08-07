#! /usr/bin/perl -w
#
# File: snp.correctString.pl
# Time-stamp: <03-Jul-2009 13:30:59 tdo>
# $Id: $
#
# Copyright (C) 2008 by Pathogene Group, Sanger Center
#icornLib::
# Author: Thomas Dan Otto
#
# Description:
#
# tdo 21.01: include read pairs of different length for SNPoMatic
#

use Data::Dumper;
use strict;
use warnings;
use lib "$ENV{ICORN_HOME}";
use icornLib;




my $string    = shift;
my $snp       = shift;
my $del       = shift;
my $ins       = shift;
my $fastq     = shift;
my $resultName= shift;

my $readPair  = shift;

my $fastqLib1 = shift;
my $lengthLib1= shift;
my $insertLib1= shift;

my $fastqLib2 = shift;
my $lengthLib2= shift;
my $insertLib2= shift;

##
my $shiftFile = shift;   # old onae
my $version   = shift;
#version=1 old, no correction

#version=2 correction with SNPoMatic, but just accept one lib with one
#length
my $color     = shift;

my $minQual=shift;
if (!defined($minQual)) {
  $minQual=60;
}
## Variable for version 2, including the shift and control with SNPoMatic
my $pileupRef = $string.".ref.fa";
my $pileupQry = $string.".qry.fa";




#version=3 accepct 2libs 190 fragement size
if (!defined($version)) {
  $version=3
}

print "version is $version Readpair = $readPair\n";


my ($ref_snp,$ref_coverage)   = icornLib::getSNP($snp,$minQual);

my $ref_del;
($ref_del,$ref_coverage)      = icornLib::getDel($del,$ref_coverage);

my $ref_fastq                 = icornLib::loadfastq($fastq);
my $ref_ins;
($ref_ins,$ref_coverage)      = icornLib::getIns($ins,$ref_fastq,$ref_coverage);

#get reference fasta
my $ref_str                   = icornLib::getFasta($string);


if ($version==1) {
  $ref_str      = icornLib::correctSequence($ref_str,
												 $ref_snp,
												 $ref_del,
												 $ref_ins);
  icornLib::writeString($ref_str,$resultName);

}
elsif ($version == 2) {
  # Do the irstcorrection
  $ref_str      = icornLib::correctSequence($ref_str,
												 $ref_snp,
												 $ref_del,
												 $ref_ins);

  icornLib::writeString($ref_str,"$resultName.tmp");

  # Generate SNPoMatic pileups of of both sequences
  icornLib::callSNPoMatic($string,$pileupRef,$fastq);   # in itereation might
                                              # not be needed, but
                                              # more robuts
  icornLib::callSNPoMatic("$resultName.tmp",$pileupQry,$fastq);
  
  #reload the reference string, as it was changed:
  my $ref_str2   = icornLib::getFasta($string);
  
  my $ref_pileupRef = icornLib::getPileupSNPoMatic($pileupRef);
  
  my $ref_pileupQry = icornLib::getPileupSNPoMatic($pileupQry);
  
  ## later the shift, as just for stats
  my $ref_sizeChr   = icornLib::getSizeChr($ref_pileupRef);
  
  my $ref_shift     = icornLib::getShift($ref_del,$ref_ins,$ref_sizeChr);
  
  $ref_str2      = icornLib::correctSequence($ref_str2,
												 $ref_snp,
												 $ref_del,
												 $ref_ins,
												 $ref_pileupRef,$ref_pileupQry,
												 $ref_shift);

  icornLib::writeString($ref_str2,$resultName);

}
#### read pair with 2 libs
elsif ($version == 3) {
  # Do the irstcorrection
  my ($ref_stats,$ref_shiftCoordinates);

  ($ref_str,$ref_stats,$ref_shiftCoordinates)
	= icornLib::correctSequence($ref_str,
									 $ref_snp,
									 $ref_del,
									 $ref_ins);
  
  icornLib::writeString($ref_str,"$resultName.tmp");

  if ($readPair && defined($lengthLib2)) {
	# Generate SNPoMatic pileups of of both sequences

	if ( -f $pileupRef ) {
	  warn "$pileupRef exists\n";
	  
	}
	else {
	  icornLib::callSNPoMatic2Paired($string,$pileupRef,
										  $fastqLib1,$lengthLib1,$insertLib1,
										  $fastqLib2,$lengthLib2,$insertLib2
										  
										 );   # in itereation might
	  # not be needed, but
	  # more robuts
	}
	
	if ( -f $pileupQry ) {
	  warn "$pileupQry exists\n";
	}
	else {
	  icornLib::callSNPoMatic2Paired("$resultName.tmp",$pileupQry,
										  $fastqLib1,$lengthLib1,$insertLib1,
										  $fastqLib2,$lengthLib2,$insertLib2);
	}
	
  }
  elsif (! $readPair) {
	icornLib::callSNPoMatic($string,$pileupRef,$fastq);   # in itereation might
	# not be needed, but
	# more robuts
	icornLib::callSNPoMatic("$resultName.tmp",$pileupQry,$fastq);
  }
  
  #reload the reference string, as it was changed:
  my $ref_str2   = icornLib::getFasta($string);
  
  my $ref_pileupRef = icornLib::getPileupSNPoMatic($pileupRef);
  
  my $ref_pileupQry = icornLib::getPileupSNPoMatic($pileupQry);
  
  ## later the shift, as just for stats
  my $ref_sizeChr   = icornLib::getSizeChr($ref_pileupRef);

    # OldShift is will be needed to map the errors to the first version
  # of the genome.
  my $ref_oldShift;
  
  if (-f "$shiftFile") {
	print "Found old shift fiel file!\n";
	
	$ref_oldShift  = icornLib::getOldShift($shiftFile);
  }
  

  #REf_shift will be need to relationate the coverage of SNPoMatic
  # between two iterations Attention $ref_shift holds all the possible
  # indels, not just the confirmed, as will be hold in $ref_shiftCoordinates
  my $ref_shift     = icornLib::getShift($ref_del,$ref_ins,$ref_sizeChr);

  my $ref_changes;
  
  ($ref_str2,$ref_stats,$ref_shiftCoordinates)
#  ($ref_str2, $ref_changes)
	            = icornLib::correctSequence($ref_str2,
												 $ref_snp,
												 $ref_del,
												 $ref_ins,
												 $ref_pileupRef,$ref_pileupQry,
												 $ref_coverage,$ref_shift);

  icornLib::writeString($ref_str2,$resultName);
  
# double check  $ref_str2   = icornLib::getFasta($string);

  icornLib::writeStats($ref_stats,$resultName,$ref_coverage,$ref_str2,
							$ref_oldShift,$color,
							$ref_pileupRef,$ref_pileupQry);
  
  my $ref_newShift = icornLib::calculateNewShift($ref_oldShift,$ref_shiftCoordinates);

  if (-f "$shiftFile") {
	print "save new\n";
	
	icornLib::saveShift($resultName,$ref_newShift);
  }
  else {
	icornLib::saveShift($resultName,$ref_shift,$ref_shiftCoordinates,$ref_sizeChr);	
  }
  

#  icornLib::writeStats($ref_stats,$resultName,$ref_coverage,$ref_str2,$ref_oldShift,$color);
  

}
#### get the stats back: This include an shifting error, how the show
# are An gff file, pointing the type of change and the coverage of the
# plot there. We can get that from the SNP and 
elsif ($version == 4) {

  #pileups do already exists.
  my $ref_pileupRef = icornLib::getPileupSNPoMatic($pileupRef);
  
  my $ref_pileupQry = icornLib::getPileupSNPoMatic($pileupQry);
  
  ## later the shift, as just for stats
  my $ref_sizeChr   = icornLib::getSizeChr($ref_pileupRef);

  # OldShift is will be needed to map the errors to the first version
  # of the genome.
  my $ref_oldShift;
  
  if (-f "$shiftFile") {
	print "Find old file!\n";
	
	$ref_oldShift  = icornLib::getOldShift($shiftFile);
  }
  

  #REf_shift will be need to relationate the coverage of SNPoMatic
  # between two iterations Attention $ref_shift holds all the possible
  # indels, not just the confirmed, as will be hold in $ref_shiftCoordinates
  my $ref_shift     = icornLib::getShift($ref_del,$ref_ins,$ref_sizeChr);

  my ($ref_stats,$ref_shiftCoordinates);
  
  ($ref_str,$ref_stats,$ref_shiftCoordinates)
	= icornLib::correctSequence($ref_str,
									 $ref_snp,
									 $ref_del,
									 $ref_ins,
									 $ref_pileupRef,$ref_pileupQry,
									 $ref_coverage,$ref_shift);
  
  my $ref_str2   = icornLib::getFasta($string);

  icornLib::writeStats($ref_stats,$resultName,$ref_coverage,$ref_str2,
							$ref_oldShift,$color,
							$ref_pileupRef,$ref_pileupQry);

  my $ref_newShift = icornLib::calculateNewShift($ref_oldShift,$ref_shiftCoordinates);
  if (-f "$shiftFile") {
	print "save new\n";
	
	icornLib::saveShift($resultName,$ref_newShift);
  }
  else {
	icornLib::saveShift($resultName,$ref_shift,$ref_shiftCoordinates,$ref_sizeChr);	
  }
}

exit 0;
