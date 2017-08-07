package icornLib;

  
# Copyright (C) 2008 by Pathogene Group, Sanger Center
#
# Author: Thomas Dan Otto
#
# Description:
#

# 18.04. ins was not the correct position. Include most frequente base
# 27.04. delete position now taken from the most frequent position
# 18.06. adapted program for standalone; bug with SN


use Data::Dumper;
use strict;

my $SNPOMATIC_PATH    = $ENV{SNPOMATIC_HOME};
my $SNPOMATIC_PROGRAM = "findknownsnps";


my $MIN_INDEL_COVER   = 10;
my $MIN_INDEL_PERCENT = (3/10);


my $debug             = 0;


#### Function that will correct the string itself.
# The $ref_str{chr}[pos] will be changed with the SNP and indels
#
#
#/
sub correctSequence{
  my $ref_str   = shift;
  my $ref_snp   = shift;
  my $ref_del   = shift;
  my $ref_ins   = shift;

  my $Version2  = 0;
  
  # for version 2
  my $ref_pileupRef = shift;
  my $ref_pileupQry = shift;
  my $ref_coverage  = shift;
  my $ref_baseDiff   = shift;
  
  if (defined($ref_baseDiff)) {
	$Version2       = 1;	
  }

  my %shiftCoordinates;
  
  
  my %Statistics;
  
  # Goes through all chromosmes and check if there exists SNP, del or
  # ins. 
  #
  # if $ref_pileup is set, which means, comparison after before is
  # available. The check is encoded in the function: isCorrectionGood
  
  foreach my $chr (keys %$ref_str) {
	
	print "working on $chr\n", if $debug;
	
	#SNP
	foreach my $pos (sort {$a <=> $b} keys %{$$ref_snp{$chr}}) {
	  # if it is not Version xor Corection is not godd do nothing
	  
	  if (!($Version2 && !(isCorrectionGood($ref_pileupRef,
											$ref_pileupQry,$ref_baseDiff,$pos,$chr)))) {
		
		if ( uc($$ref_str{$chr}[($pos-1)]) ne uc($$ref_snp{$chr}{$pos})
			 && !(defined($$ref_ins{$chr}{$pos}))
			 && !(defined($$ref_del{$chr}{$pos}))
			 
		   ) {
		  $Statistics{SNP}{$chr}{($pos-1)} = $$ref_str{$chr}[($pos-1)];
		  $$ref_str{$chr}[($pos-1)] = $$ref_snp{$chr}{$pos};
		  print "Changing pos $pos on $chr\n", if $debug;
		}
		elsif (defined($$ref_ins{$chr}{$pos})
			   || defined($$ref_del{$chr}{$pos})) {
		  
		}
		else {
		  if (defined($$ref_ins{$chr}{$pos})
			  || defined($$ref_del{$chr}{$pos})) {
			
		  }
		  else {
			$Statistics{HETERO}{$chr}{($pos-1)} = $$ref_snp{$chr}{$pos};
		  }
		}
		
	  }
	  ### log how many SNP were rejected.
	  else {
		$Statistics{"Rej.SNP"}{$chr}{($pos-1)}=$$ref_str{$chr}[($pos-1)];
	  }
	}
#		print Dumper %{$$ref_del{$chr}};
  
	# DELETION
	foreach my $pos (sort {$a <=> $b} keys %{$$ref_del{$chr}}) {
	  if (!($Version2 && !(isCorrectionGood($ref_pileupRef,
											$ref_pileupQry,$ref_baseDiff,$pos,$chr)))) {
		if (! defined($$ref_ins{$chr}{$pos})) {

		  my $size = $$ref_del{$chr}{$pos};
	# not the size, just bases $Statistics{DEL}{$chr}{($pos-1)} = $size;
		  $shiftCoordinates{$chr}{($pos-2)} = ((-1)*$size);
		  
		  
		  # The deletion will be from pos-1, with the siye $size:
		  my $posToDelete = ($pos-1);
		  
		  print "Changing pos $pos on $chr\n", if $debug;
		  
		  while ($size>0) {
			$Statistics{DEL}{$chr}{($pos-1)}.= $$ref_str{$chr}[$posToDelete];
			$$ref_str{$chr}[$posToDelete] ='';
			
			
			$size--;
			$posToDelete++;
		  }
		}
		else {
		  print "Changing pos $pos on $chr\n", if $debug;
		  print "Ignored deletion $chr $pos \n";
		}
		
	  }
	  ### log how many SNP were rejected.
	  else {
		$Statistics{"Rej.DEL"}{$chr}{($pos-1)}=$$ref_del{$chr}{$pos};
	  }

	}
	
		print "INSERTION\n", if $debug;
	#INSERTION
#	print Dumper %{$$ref_ins{$chr}};
  
	foreach my $pos (sort {$a <=> $b} keys %{$$ref_ins{$chr}}) {
	  
	  if (!($Version2 && !(isCorrectionGood($ref_pileupRef,
											$ref_pileupQry,$ref_baseDiff,$pos,$chr)))) {

		
		$Statistics{INS}{$chr}{($pos-1)}.=$$ref_str{$chr}[($pos-1)];
		$$ref_str{$chr}[($pos-1)] =$$ref_ins{$chr}{$pos}.$$ref_str{$chr}[($pos-1)];
		$shiftCoordinates{$chr}{($pos-2)}=(length($$ref_ins{$chr}{$pos}));
		
		print "Inserted changed:  Changing pos $pos on $chr\n", if $debug;

	  }
	  ### log how many SNP were rejected.
	  else {
		print "Rejected Inserted changed:  Changing pos $pos on $chr\n", if $debug;
		$Statistics{"Rej.INS"}{$chr}{($pos-1)}.=$$ref_ins{$chr}{$pos};
	  }
	  
	}
  }
  
  return ($ref_str,\%Statistics,\%shiftCoordinates);
}



### Will check is Ref is smaller than Qry, and return 1.  If both 0
# will also return 1, as the position could be in process of
# correction TODO Actually, if the coverage is lower, it would be
# interesting to see, what the so far covering reads are doing. If
# they are repeats or 
sub isCorrectionGood{
  my $ref_pileupRef = shift;
  my $ref_pileupQry = shift;
  my $ref_baseDiff   = shift;
  my $pos = shift;
  my $chr = shift;
  

  my $return =0;
  
  my $covRef=$$ref_pileupRef{$chr}[$pos];
  my $posQry=($pos+$$ref_baseDiff{$chr}[$pos]);
	  
  my $covQry=$$ref_pileupQry{$chr}[$posQry];

  print "$pos : $covRef :  $posQry : $covQry : $$ref_baseDiff{$chr}[$pos] \n", if ($debug);

  print "$covRef > $covQry\n", if $debug;
  
  if ($covRef > $covQry) {
	print "BAD\n", if $debug;
	
	$return= 0;
  }
  else {
	$return = 1;
  }
  return $return;
}


#### Function that will write the ref_str hash to a multifaste file.
# INPUT: $ref_str{chr}[pos] and $resutName.
#/
sub writeString{
  my $ref_str   = shift;
  my $resultName= shift;
  
  
  open (F,"> $resultName") or die "Problem to save $resultName: $!\n";
  foreach my $chr (sort keys %$ref_str) {
	print F ">$chr\n";
	for (my $i=100; $i<scalar(@{$$ref_str{$chr}});$i+=100) {
	  $$ref_str{$chr}[$i].="\n";
	}
	print F join ("",@{$$ref_str{$chr}}),"\n";
  }
  
  close(F);
}

sub getOldShift{
  my $filename = shift;


  open (F, $filename ) or die "Couldn't find fie $filename \n";

  my @F =<F>;
  close(F);
  
  my $name;
  my %shift;
  
  foreach (@F) {
	if (/>(\S+)/) {
	  $name = $1;
	}
	else {
	  chomp;
	  push @{ $shift{$name}}, $_;
	  
	}
  }
  close(F);
  
  return (\%shift);
}

sub calculateNewShift{
  my $ref_oldShift    = shift;
  my $ref_shiftCoordinates = shift;

  my %ref_newShift;
  
  foreach my $chr (keys %$ref_oldShift) {
	my $shiftValueNew=0;
	for (my $i=0; $i < scalar(@{$$ref_oldShift{$chr}});$i++) {
	  if (defined($$ref_shiftCoordinates{$chr}{$i})) {
		$shiftValueNew+=$$ref_shiftCoordinates{$chr}{$i};
	
		
	  }
	  push @{ $ref_newShift{$chr} }, ($$ref_oldShift{$chr}[$i]+$shiftValueNew);
	}
  }
  return (\%ref_newShift)
}

sub saveShift{
  my $resultName   = shift;
  my $ref_newShift = shift;
  my $ref_shiftCoordinates = shift;
  my $ref_sizeChr          = shift;
  
  
  my $res;
  # if ture, them just save the shift
  if (!defined($ref_shiftCoordinates)) {
  
	foreach my $chr (keys %$ref_newShift) {
	  $res .= ">".$chr."\n";
	  foreach (@{$$ref_newShift{$chr}}) {
		$res .= "$_\n";
	  }
	}
  }
  # else , construct the shift
  else {
	
	foreach my $chr (keys %$ref_sizeChr) {
	  $res .= ">".$chr."\n";
	  my $shiftValueNew=0;
	  
	  for (1..($$ref_sizeChr{$chr})) {
		if (defined($$ref_shiftCoordinates{$chr}{$_})) {
		  $shiftValueNew+=$$ref_shiftCoordinates{$chr}{$_};
		 

		}
		$res .= ($shiftValueNew);
		$res .= "\n";
	  }
	}
  }
  
  open (F, "> $resultName.shift.txt") or die "Couldn't open file $resultName.shift.txt";
  print F $res;
  close(F);
}

#### Function that will write the ref_str hash to a multifaste file.
# INPUT: $ref_str{chr}[pos] and $resutName.
#/
sub writeStats{
  my $ref_stats   = shift;
  my $resultName= shift;
  my $ref_coverage = shift;
  my $ref_str      = shift;
  my $ref_oldShift = shift;
  my $color     = shift;
  my $ref_pileupRef = shift;
  my $ref_pileupQry = shift;
  
  $color=($color%15);
  
  
  my $resPerBase;
  my $resGeneral;
  my $ref_resGFF;

  my @show=('SNP','INS','DEL','HETERO','Rej.SNP','Rej.INS','Rej.DEL');
  
#  foreach my $what (sort keys %$ref_stats) {
  foreach my $what (@show) {
	
	$resPerBase .= "$what\n";
	$resGeneral .= "$what\n";
	foreach my $chr (sort keys %$ref_str) { # ref_str has all replicons
	  $resPerBase .= "$chr\n";
	  $resGeneral .= "$chr\t";
	  my $count=0;
	  foreach my $pos (sort {$a <=> $b} keys %{$$ref_stats{$what}{$chr}}) {


		my $posIter = ($pos);
		my $posRef  = ($pos+1);
		if (defined($$ref_oldShift{$chr}[$pos])) {
		  $posRef = ($pos+1-$$ref_oldShift{$chr}[$pos]);
		  
		}
		
		$resPerBase .= "$chr\t$pos\t(".($posRef-1).")\t$what\t";
		if (!defined($$ref_stats{$what}{$chr}{$pos})) {
		  $resPerBase .="\n";
		  print "Problem with $chr $what $pos \n";
		  die;
		  
		}
		else {		
		  $resPerBase .=  $$ref_stats{$what}{$chr}{$pos}." -> $$ref_str{$chr}[$pos] \n";;
		}
		### As positinos are already ok! oiin array start 0 get a -1
		### to each pos

		

		$ref_resGFF = write2GFF($ref_stats,$ref_resGFF,$chr,$what,$pos,$posRef,$posIter,$ref_coverage,$$ref_str{$chr}[$pos],$color,$ref_pileupRef,$ref_pileupQry);
		
		
		$count++;
	  }
	  $resGeneral .= "$count\n";
	}

	foreach my $chr (keys %$ref_resGFF) {
	  open (F,"> $resultName.$chr.gff") or
	  die "Problem to save $resultName.General.stats: $!\n";
		print F $$ref_resGFF{$chr};
		close(F);
	  
	}

	open (F,"> $resultName.General.stats") or
	  die "Problem to save $resultName.General.stats: $!\n";
	print F $resGeneral;
	close(F);
	
	open (F,"> $resultName.PerBase.stats") or
	  die "Problem to save $resultName.PerBase.stats: $!\n";
	print F $resPerBase;
	close(F);
  }
}

sub write2GFF{
  my $ref_stats = shift;
  my $ref_resGFF= shift;
  my $chr       = shift;
  my $what      = shift;
  my $pos       = shift;
  my $posRef    = shift;
  
  my $posNow    = shift;
  
  my $ref_coverage = shift;
  my $basenow   = shift;
  
  my $color     = shift;

  my $ref_pileupRef = shift;
  my $ref_pileupQry = shift;
  ### 

#  $$ref_resGFF{$chr} .= "unknown\tCor\tCDS\t".($posRef)."\t".($posRef)."\t0\t+\t.\tnote \"$what; Coverage: ".$$ref_coverage{$chr}{($pos+1)}."  ".$$ref_stats{$what}{$chr}{$pos}." Base in iteration before: $basebefore; Pos in this iteration $posNow\" ; color=$color label=\"$what\"\n";
  chomp($basenow);
  
  $$ref_resGFF{$chr} .= "unknown\tCor\tCDS\t".($posRef)."\t".($posRef)."\t0\t+\t.\tcolor=$color;label=\"$what\";note=$what:+Coverage:".$$ref_coverage{$chr}{($pos+1)};
  $$ref_resGFF{$chr} .= "++New+BASE:+$basenow++before:+".$$ref_stats{$what}{$chr}{$pos}."++Pos+in+this+iteration+$posNow+++CoveragePerfect (before):+";
  $$ref_resGFF{$chr} .= $$ref_pileupRef{$chr}[($pos+1)]."+++(now):+".$$ref_pileupQry{$chr}[($pos+1)]."+++Iteration+$color\"\n";
  
  # print "unknown\tCor\tCDS\t$pos\t$pos\t0\t+\t.\tnote \"$what Coverage
  #".$$ref_coverage{$chr}{$pos}."  ".$$ref_stats{$what}{$chr}{$pos}." 
  #Base in iteration before: $basebefore; Pos in this iteration
  #$posNow\" ; label=\"$what\"\n"; print
  #"unknown\tCor\tCDS\t$pos\t$pos\t0\t+\t.\tnote \"$what: Coverage
  #".$$ref_coverage{$chr}{$pos}." ".$$ref_stats{$what}{$chr}{$pos}." 
  #\" ; c label=\"$what\"\n";
  return $ref_resGFF;
}

sub getDel{
  my $name=shift;
  my $ref_coverage = shift;
  
  
  my @ar2;
  open (F, "$name") or die "probelms $name\n";
  
  
  @ar2=<F>;
  close(F);
  
  my %h;
  # print Dumper @ar2;
  
  
  for (my $i=3; $i<(scalar(@ar2)-5);$i++) {
	chomp($ar2[$i]);
	my @ar=split(/\s+/,$ar2[$i]);
	#	print $ar2[$i];
	
	my $count=0;
	my $size=0;
	my %hash_position;
	
	while (defined($ar[0]) && $ar[0] eq 'Deletion:') {
	  $count++;
	  $size+=$ar[4];
	  $i++;
	  
	  chomp($ar2[$i]);
	  $hash_position{$ar[3]}++;

	  #	  if (defined($$ref_ar[$i]) and $$ref_ar[$i] ne '') {
	  @ar=split(/\s+/,$ar2[$i]);
	 
	}

	@ar=split(/\s+/,$ar2[($i-1)]);

		   
#	print $ar2[($i-2)]."\n";
#	if (!defined($count) or !defined($ar[10]) or !defined($MIN_INDEL_COVER)) {
#	  print $ar2[($i-1)]."\n";
#	  print " $MIN_INDEL_COVER \n";
#	  print " $count :: $ar[10] \n";
#	  
#	  
#	}
	
	if ($count >= $MIN_INDEL_COVER &&
		$count > ($MIN_INDEL_PERCENT*$ar[10])) {
	  
	  my $chr      = $ar[2];
	  my $coverage = $ar[10];

	  my @test;
	  foreach  (reverse sort {$hash_position{$a} <=> $hash_position{$b}} keys %hash_position){
		push @test, $_;
		
	  }

	  my $pos      = $test[0];
	  
	  
	  $$ref_coverage{$chr}{$pos}+=$coverage;
	  #		print $ar2[($i-2)]."\n",if $debug;	
	  # $ar[4] holds the size of the deletion
	  $h{$chr}{$pos}=$ar[4];
	}
  }
  
  return (\%h,$ref_coverage);
}


sub getSNP{
  my $name = shift;
  my $minQual = shift;
  
  
  open (F, "$name") or die "probelms $name\n";

  #Will hold the coverage of each SNP and indel
  my $ref_coverage;
  
  my %h;
  
  my @ar=<F>;
  close(F);
  
  my @res;
  
  ## plot of SNP
  for (my $i=7;$i<scalar(@ar);$i++) {
	chomp($ar[$i]);
	my @ab =split(/\s+/,$ar[$i]);
	#SNP_hez: Tbg972_01; 10 3         4       A C/A 1      2      1      0      0      0      0   0   0   0
	if (defined($ab[9]) and $ab[2] >= $minQual) {
	  $ab[1]=~ s/;//g;

	  
	  my $chr      = $ab[1];
	  my $coverage = $ab[4];
	  my $pos      = $ab[3];
	  
	  
	  $$ref_coverage{$chr}{$pos}+=$coverage;

	  my %small;
	  $small{($ab[7]-$ab[13])}='A';
	  $small{($ab[8]-$ab[14])}='C';
	  $small{($ab[9]-$ab[15])}='G';
	  $small{($ab[10]-$ab[16])}='T';
	  my @test;
	  
	  foreach  (reverse sort {$a <=> $b} keys %small){
		push @test, $_;
	  }

	  
	  $h{$ab[1]}{$ab[3]}=$small{$test[0]};
	}
  }
  return (\%h,$ref_coverage);
}




sub getFasta{
  open F, shift or die "prob couldn't find fiil  \n";

  my %h;
  my @ar;
  my $name;
  
  while (<F>) {
        chomp;
        if (/^>(\S+)/){
		  $name=$1;
		}
		else {
		  chomp;
		  my @a=split(//);
		  foreach (@a) {
			push @{$h{$name}}, $_;
		  }
		}
	  }
  close(F);
  
  return \%h;
}

sub getIns{
  my $name=shift;
  my $ref_fastq=shift;
  my $ref_coverage = shift;
  
  
  open (F, "$name") or die "couldn't open file $name \n";
  
  my @File=<F>;
  my $ref_File=\@File;
  close(F);
  
  my @ar;
  
  my %h;
  
  
  for (my $i=3; $i<(scalar(@$ref_File)-5);$i++) {
	chomp($$ref_File[$i]);
	my @ar=split(/\s+/,$$ref_File[$i]);
	
	my $count=0;
	my $size=0;
	my $bases='';
	my $base;
	my %hash_bases;
	my %hash_pos;
	
	
	while (defined($ar[0]) && $ar[0] eq 'Insertion:') {
	  $count++;
	  $size+=$ar[4];
	  $i++;
	  chomp($$ref_File[$i]);

	  # Attation: 0 - formward and 1 reverse... SSAHA pileup stuff
	  if ($ar[9]) {
		$base = revcomp(substr($$ref_fastq{$ar[7]},($ar[8]-1),$ar[4]));
	  }
	  else {
		$base =substr($$ref_fastq{$ar[7]},($ar[8]-1),$ar[4]);
	  }
	  
	  $hash_bases{$base}++;
	  $hash_pos{$ar[3]}++;
	  
	  @ar=split(/\s+/,$$ref_File[$i]);
	  
	}
	@ar=split(/\s+/,$$ref_File[($i-5)]);
	
	
	
	if ($count >= $MIN_INDEL_COVER &&
		$count > ($MIN_INDEL_PERCENT*$ar[10])) {
	
	  my $chr      = $ar[2];
	  my $coverage = $ar[10];
	### this changed to get best postion  my $pos      = $ar[3];
	  

	  my @getPos;
	   foreach  (reverse sort {$hash_pos{$a} <=> $hash_pos{$b}} keys %hash_pos){
		push @getPos, $_;
	  } 
	  my $pos = $getPos[0];
	  $$ref_coverage{$chr}{$pos}+=$coverage;

#	  print "$ar[2] Pos is now $pos \n";
	  
	  
	  my @test;
	  foreach  (reverse sort {$hash_bases{$a} <=> $hash_bases{$b}} keys %hash_bases){
		push @test, $_;
	  }

	  if ($ar[9]) { ## reverse
		$base = revcomp(substr($$ref_fastq{$ar[7]},($ar[8]-1),$ar[4]));
	  }
	  else {
		$base = substr($$ref_fastq{$ar[7]},($ar[8]-1),$ar[4]);
	  }
	  
	  $h{$ar[2]}{$pos}=$test[0];
	}
  }
 
  return (\%h,$ref_coverage)
}

sub countIns{
  my $name=shift;

  open (F, "$name") or die "couldn't open file $name \n";

  my @array=<F>;
  my $ref_File=\@array;
  close(F);
  
  my @ar;
  
  my %h;
  
  
  for (my $i=3; $i<(scalar(@$ref_File)-5);$i++) {
	chomp($$ref_File[$i]);
	my @ar=split(/\s+/,$$ref_File[$i]);
	
	my $count=0;
	my $size=0;
	my $bases='';
	my $base;
	
	while (defined($ar[0]) && $ar[0] eq 'Insertion:') {
	  $count++;
	  $size+=$ar[4];
	  $i++;
	  chomp($$ref_File[$i]);
	  @ar=split(/\s+/,$$ref_File[$i]);
	}
	@ar=split(/\s+/,$$ref_File[($i-1)]);
	
	if ($count > (0.4*$ar[10])) {
	  $h{$ar[2]}{$ar[3]}="";
	  
	}
	
  }
  return (\%h)
}

sub loadfastq{
  my $name=shift;
  
  my %h;
  keys( %h )=1000000;
  
  open (F,$name) or die "problem to load fastq $name: $! \n";

  my $n;
  
  while (<F>) {
	chomp;
	($n) = $_ =~ /^@(\S+)/;
	$_=<F>;
	chomp;
	$h{$n}=$_;
	$_=<F>;	
	$_=<F>;
  }
  close(F);
  
  return \%h;
}


sub revcomp{
  my $rev=reverse(shift);
  
  $rev =~ tr/ATGC/TACG/;
  return $rev;
  

}

#### Functino to load the pileup in %h{chr}[pos]=coverage
# format of line is:
# MAL1    C       7       C       20      c cc      c  c    C  cCC    cC   C   C CC   c     c      c   c     C    ...
#
  # Input: pileupfile
  ##/
sub getPileupSNPoMatic{
  my $name=shift;

  open (F, "$name") or die "couldn't open file $name \n";

#  my @array=<F>;

  my %h;

  # check if maybe the wrong format is used:
#  my @ar=split(/\t/,$array[10]);
#  if ((defined($ar[5]) || defined($ar[4]))) {
#	print "Sorry, wrong SNPoMatic format in $array[10]\n";
#	print "please do a awk '{ print $1\"\t\"$3\"\t\"$5  }'  3D7.genome.1.pileup... \n";
#	die "Exits in getPileupSNPoMatic\n"
#  }
  
#  foreach (@array) {
  while (<F>) {
  chomp;
	
	my @ar=split(/\t/);
	if (!defined($ar[2])) {
	  $ar[2]=0;
	  
	}
	$h{$ar[0]}[$ar[1]]=$ar[2];
  }

  print "Have pileup done $name.\n",if $debug;
  close(F);
  
  return \%h;
}

#### Functnio to get the length of each chromosom
# INPUT: $ref_pileup of the reference
# OUTPUT: %h{chr}=Length
##/
sub getSizeChr{
  my $ref_pileup = shift;
  my %h;
  
  foreach my $chr (keys %$ref_pileup) {
	$h{$chr}=(scalar(@{$$ref_pileup{$chr}})-1);	
  }
  return (\%h)
}
  
### Function to get the difference of each position between two iterations
# Through the indels this can be detected easely.
# Each shift (indel) is relative to the reference. Ins is +1 a Del -1
#
# INPUT: $ref_ins{chr}{pos} and $ref_del{chr}{del} $ref_sizeChr{chr}
#
# OUTPUT: %h{chr}[pos]=shift (integer)
#
#/
sub getShift{ 
  my $ref_del = shift;
  my $ref_ins = shift;
  my $ref_sizeChr = shift;
  

  my %h;
  
  foreach my $chr (sort keys %$ref_sizeChr) {
	my $diff=0;

	print "$chr $$ref_sizeChr{$chr} \n", if $debug;
	
	# run through each chromooms
	foreach my $pos (0..$$ref_sizeChr{$chr}) {

	  # if finds ins, $diff+= size of it
	  if (defined($$ref_ins{$chr}{$pos})) {
		$diff+= length($$ref_ins{$chr}{$pos});
		print "Found insert $diff\n", if $debug;		
	  }

	  # if finds del, $diff-= size of it
	  if (defined($$ref_del{$chr}{$pos})) {
		$diff-= $$ref_del{$chr}{$pos};
		print "Found deletion $diff\n", if $debug;
		
	  }

	  #Write the diference
	  $h{$chr}[$pos]=$diff;
	  
	}
  }
  print "Shifted done.\n",if $debug;
  
  return \%h;
}


### Will call findknownsnp, defeind in $SNPOMATIC_PATH and
# $SNPOMATIC_PROGRAM Generated will be a pileup file. As this file is
# to big, it is afterwards awk:
#
#
# for futher processing. The old pilup will be deelted.
sub callSNPoMatic{
  my $ref = shift;
  my $pileup = shift;
  my $fastq = shift;

  my $tmp_pileup = "$pileup.tmp";
  
  my $Call ="$SNPOMATIC_PATH/$SNPOMATIC_PROGRAM --genome=$ref --fastq=$fastq --pileup=$tmp_pileup";
  !system("$Call") or die "SNPoMatic did not run ok $Call \n";

  
  $Call="awk '{ print \$1\"\\t\"\$3\"\\t\"\$5 }'  $tmp_pileup  > $pileup ";
  
  !system("$Call") or die "SNPoMatic did not run ok $Call \n";

  unlink($tmp_pileup);
  
}
### Will call findknownsnp, defeind in $SNPOMATIC_PATH and
# $SNPOMATIC_PROGRAM Generated will be a pileup file. As this file is
# to big, it is afterwards awk:
#
# The difference to the previous function is the possibilty to go for
# pairs with different length
#
# for futher processing. The old pilup will be deelted.
sub callSNPoMatic2Paired{
  my $ref        = shift;
  my $pileup     = shift;
  my $fastqLib1  = shift;
  my $lengthLib1 = shift;
  my $insertLib1 = shift;
  
  my $fastqLib2  = shift;
  my $lengthLib2 = shift;
  my $insertLib2 = shift;
  

  callSubSNPoMaticPaired($ref,"$pileup.Lib1",$fastqLib1,$lengthLib1,$insertLib1);
  callSubSNPoMaticPaired($ref,"$pileup.Lib2",$fastqLib2,$lengthLib2,$insertLib2);

  my $Call="paste -d: $pileup.Lib1 $pileup.Lib2  > $pileup.tmp";

  
  !system("$Call") or die "awk went wrong did not run ok $Call \n";

  ### the $pileup.tmp will hold two 3 line pileups. We want the first
  ### three and the last...
  $Call="awk -F':' '{ print \$1\"\\t\"\$2\"\\t\"\$3+\$6 }'  $pileup.tmp  > $pileup";
  
  !system("$Call") or die "SNPoMatic did not run ok $Call \n";

  unlink("$pileup.Lib1");
  unlink("$pileup.Lib2");
  unlink("$pileup.tmp");
  
}

sub callSubSNPoMaticPaired {
  my $ref        = shift;
  my $pileup     = shift;
  my $fastq      = shift;
  my $length     = shift;
  my $insert_size= shift;
  
  my $tmp_pileup = "$pileup.$length.tmp";
  
  # call SNP
  my $Call ="$SNPOMATIC_PATH/$SNPOMATIC_PROGRAM --genome=$ref --fastq=$fastq --pileup=$tmp_pileup --pair=$length --fragment=$insert_size --chop=2";
  !system("$Call") or die "SNPoMatic did not run ok $Call \n";
  
  $Call="awk '{ print \$1\":\"\$3\":\"\$5 }'  $tmp_pileup  > $pileup";
  
  !system("$Call") or die "awk did not run ok $Call \n";
  
  unlink($tmp_pileup);
}

1;
