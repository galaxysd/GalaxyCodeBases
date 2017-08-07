#!/bin/bash

# # In no event shall the authors of the software or GRL be responsible #
# # or liable for any loss or damage whatsoever arising in any way      #
# # directly or indirectly out of the use of this software or its       #
# # derivatives, even if advised of the possibility of such damage.     #
# # Our software can be freely distributed under the conditions set out #
# # above, and must contain this copyright notice.                      #
# #######################################################################
# *
# *  Author : Thomas D. Otto (tdo@sanger.ac.uk)
# *
# *  Copyright (C) 2009 by Genome Research Limited, All rights reserved.
# *
# **************************************************************************/
# Description:
# Acorn is an iterative correction algroithm. It can also be used to find more SNP's
# and confirm there significance in term of correct calling. First the sequences are
# adapted to the correct format. Next they are mapped with SSAHA against the REFERENCE.
# SNP's and short indels are called by ssaha pileup. The REFERENCE is corrected.
# Perfect coverage of changes is calculated, to possible revert. 
# Script finish, once no more erros can be found.
#
# Requirements:
# ssaha_pileup in path and working
# SNPomatic https://sourceforge.net/projects/snpomatic/ in path.

ICORN_HOME=/scratch/wangq4/bin/iCORN
PILEUP_HOME=/scratch/wangq4/bin/iCORN/pileup_v0.5b
SNPOMATIC_HOME=/scratch/wangq4/bin/iCORN
SsahaDoCoveragePlots=1
export ICORN_HOME PILEUP_HOME SNPOMATIC_HOME SsahaDoCoveragePlots

# will map es the script before.



# tipical cals
# nohup  ~/Bin/helper.mapTrans3D7.v2.iter.RP.sh RunPaper 3D7.sequence.2.1.4.NoN.fasta 1 1000000 1022_s_All.fastq 36 1105_s_All.fastq 76 &

# Collect input paramter
RefAll=$1
#readPair=$2;
start=$2
end=$3

fastq1=$4
fastq2=$5
insertSize=${6}
insertlib1=$7


lengthlib1=$(head -n 2 $fastq1 | tail -n 1 | perl -nle 'print length $_')
echo "readlenfth =$lengthlib1"
# STATIC VARIABLES
LSF=0;    # if set to 1, an LSF cluster is assumed.


#get just the reference.
ref=$(echo $RefAll | perl -nle '@ar=split(/\//); print ($ar[(scalar(@ar)-1)])' )
ln -s $RefAll $ref.1
 
### check if the program can be found
# snpomatic
if [ ! -f "$SNPOMATIC_HOME/findknownsnps" ] ; 
then 
	echo "Sorry the program SNPomatic couldn't be found. Please install it and make sure, the varible SNPOMATIC_HOME is set.";
	echo "But trying < ls \$SNPOMATIC_HOME/findknownsnps >"
	echo "The program should be found."
	exit;
fi

# ssaha
if [ ! -f "$PILEUP_HOME/pileup.csh" ] ; 
then 
	echo "Sorry the program SSAHA_pileup couldn't be found. Please install it and make sure, the varible PILEUP_HOME is set.";
	echo "But trying < ls \$PILEUP_HOME/pileup.csh >"
	echo "The program should be found."
	exit;
fi

# icorn
if [ ! -f "$ICORN_HOME/icornLib.pm" ] ; 
then 
	echo "Sorry the home directory for ICORN isn't set correctly. Please make sure, the varible PILEUP_HOME is set.";
	echo "But trying < ls \$ICORN_HOME/icornLib.pm >"
	echo "a file should be found."
	exit;
fi




### get the external global variable. Those do not change so often, and should be just changed if needed.

# Quality need for evidence of error.
if [ -z "$CARMA_CORRECT_QUAL" ] 
	then
	minQual=60
else
	minQual=$CARMA_CORRECT_QUAL;
fi

# split size of fastq file
if [ -z "$CARMA_SPLICE_SIZE" ] 
	then
	splitSize=900000
else
	splitSize=$CARMA_SPLICE_SIZE
fi

# Memory request if bigger genomes (>200MB 4900 is good)
if [ -z "$SsahaMemoryLimit" ] 
	then 
	memoryNeed=1900
else 
	memoryNeed=$SsahaMemoryLimit
fi



if [ -z "$fastq2" ]; then 
	readPair=0;
else 
	readPair=1;
fi


### start for itereation

echo " Calles helper.maTran $root $ref  $splitSize $lib1 $lengthlib1 $insertlib1 $lib2 $lengthlib2 $insertlib2 $start $end $insertSize $minQual";


if [ -z "$fastq1" ] && [ ! -z "$fastq2" ] ; then 
	echo "usage for single lane: $0 <reference> <iteration start> <iteratoin end 6> <fastq>";
	echo
	echo "$0 will map iteratively the sequences of the <root.mate.fastq> against the reference genome. <reference size> is the length of the reference genome. <splitSize> is the amount of fastq fisequences per chunk"
	echo 
	echo "importante:"
	echo "the fastq file must be splited in pairs, with the name <root>.mate.fastq"
	echo "The reference must be named <reference>.$start. This file must exist"
	exit
else 
	if [ -z "$insertlib1" ] ; then
		echo "usage for read pairs: $0 <reference> <iteration start> <iteratoin end 6> <fastq.F> <fastq.R> <Insert range like 100,400> <mean insert size>"
		
		echo
		echo "$0 will map iteratively the sequences of the <root.mate.fastq> against the reference genome. <reference size> is the length of the reference genome. <splitSize> is the amount of fastq fisequences per chunk"
		echo 
		echo "importante:"
		echo "the fastq file must be splited in pairs, with the name <root>.mate.fastq"
		echo "The reference must be named <reference>.$start. This file must exist"
		exit
	fi
	
fi

## set the statr and end for the amount of iteration
if [ -z "$start" ]; then
	start=1;
fi
if [ -z "$end" ]; then
	end=5
fi

if [ -z "$insertSize" ] ; then
	insertSize="70,300";
fi

echo "$start $end $otherSolexa $insertlib2"



# manipulate the fastq files...
if [ "$readPair" == "0" ] 
	then
	ln -s $fastq1 $fastq.mate.fastq
	root=$fastq1
	lib1=$fastq1
else 
	if [ -z "$fastq2" ] ; then echo "please specify second fastq file"; exit; fi; 
	echo "transform.fastq2fastqMate.pl $fastq1 $fastq2 $root"
	root=$(echo $fastq1 | sed 's/\.gz//g' | sed 's/\.F//g' | sed 's/\.fastq//g')
	echo
	echo
	$ICORN_HOME/transform.fastq2fastqMate.pl $fastq1 $fastq2 $root.mate.fastq
	ok=$?
	if [ $ok != "0" ] ; then 
		echo "Sorry something went wrong with your fastq file. Please read above!"
		exit;
	fi
	echo
	echo "transform.fastqMate2fastqSNPomatic.pl $root.mate.fastq $root.SNPomatic.fastq"
	$ICORN_HOME/transform.fastqMate2fastqSNPomatic.pl $root.mate.fastq $root.SNPomatic.fastq;
	lib1=$root.SNPomatic.fastq
	ok=$?
	if [ $ok != "0" ] ; then 
		echo "Sorry something went wrong with your fastq file. Please read above!"
		exit;
	fi
fi

mkdir tmp
cd tmp


### check if files ok
if [ ! -f "../$ref.$start" ]                # be sure the directory /mnt exists
	then
    echo "Reference $ref.$start doesn't exist"
	exit
fi


if [ ! -f "../$root.mate.fastq" ]                # be sure the directory /mnt exists
	then
    echo "mate file $root.mate.fastq doesn't exist"
    exit
fi

if [ ! -f "../$lib1" ]                # be sure the directory /mnt exists
	then
    echo "SNP-O-Matic fastq $lib1 doesn't exist"
    exit
fi


cd ..

# a debug information due to compartibility
head -n 4 $lib1 > small.fastq
lib2=small.fastq
	
	
# get an unique IF
tmpId=$$
# tranform fastq to correct

echo "Iterative call of $ref.$start for $start to $end iterations with sequence $root.mate.fastq"
echo " "

for ((i=$start;$i<=$end;i++)); do
# split the reads for ssaha in .F .R
#~tdo/bin/pileup_v0.4/ssaha_pileup/ssaha_pileup/ssaha_reads $root.fastq $root.mate.fastq
	
	mkdir $ref$i
	lfs setstripe $ref$i 0 -1 -1
	cd $ref$i
	
	if [ "$LSF" = "1" ] ; then
		## split the stuff
		num=$(perl $ICORN_HOME/splitit.pl ../$root.mate.fastq $splitSize $root)
	
		### collect by the return of the perl program, how many chunks there are
		echo "$root.mate.fastq was splitted in $num chunks"
		export num
	fi # end if LSF

	
	### call the mapping
	export root
	echo "Iteration $i of $end "
	
	nameFirstStage=$tmpId.$root.A.$i;

   ### will be just used, if not in the first iteration.
	nameSecondStage=$tmpId.$root.B.$i;
	t=$(($i-1))
	nameSecondStagePlusOne=$tmpId.$root.B.$t;
	
	if [ "$i" = "$start" ] ; then
		echo "first iteration. Please wait..."
		if [ "$LSF" = "1" ] ; then
			bsub -q normal  -J"$nameFirstStage[1-$num]" -R "select[type==X86_64 && mem > $memoryNeed] rusage[mem=$memoryNeed]" -M"$memoryNeed"000 -e out.%I.e -o out.%I.o "$ICORN_HOME/need.sh ../$ref.$i $root $readPair $insertSize"
		else # LSF = 0 
		   if [ "$readPair" = "1" ] ; then
				$ICORN_HOME/pileup.NoPile.csh -rtype solexa -insertSizeRange $insertSize ../$root.mate.fastq ../$ref.$i $root &> tmp.$i.out
			else 
    			$ICORN_HOME/pileup.NoPile.csh -rtype solexa $extra -paired 0 ../$root.mate.fastq ../$ref.$i $root &> tmp.$i.out 
			fi # end readPair
		fi # end if LSF



	else # $i > $start
		echo "Starting iteration ".$(($i-$start+1))
		
# t number of last iteration
		t=$(($i-1))
		
		if [ "$LSF" = "1" ] ; then
			bsub -q long -w"done($nameSecondStagePlusOne)"  -J"$nameFirstStage[1-$num]" -R "select[type==X86_64 && mem > $memoryNeed] rusage[mem=$memoryNeed]" -M"$memoryNeed"000 -e out.%I.e -o out.%I.o "$ICORN_HOME/need.sh ../$ref.$i $root $readPair  $insertSize"
		else #LSF = 0
			if [ "$readPair" = "1" ] ; then
				$ICORN_HOME/pileup.NoPile.csh -rtype solexa -insertSizeRange $insertSize ../$root.mate.fastq ../$ref.$i $root  &> tmp.$i.out
			else 
   				$ICORN_HOME/pileup.NoPile.csh -rtype solexa $extra -paired 0 ../$root.mate.fastq ../$ref.$i $root &> tmp.$i.out
			fi # end readPair
		fi # end if LSF =1
	fi # end if $i = $start
	
	
   ### new number of the genome.
	j=$(($i+1))
	
	if [ "$LSF" = "1" ] ; then
		bsub -q normal -w"done($nameFirstStage)" -J"$nameSecondStage" -R "select[type==X86_64 && mem > 18000] rusage[mem=18000]" -M30000000 -e out.$root.$i.e -o out.$root.$i.o "$ICORN_HOME/pileup.joinData.sh $root ../$root.mate.fastq ../$ref.$i 0 $genomesize $readPair $insertSize && $ICORN_HOME/helper2.sh $root && perl $ICORN_HOME/icorn.Correct.pl ../$ref.$i $root.snp  $root.del $root.ins Ins.Reads.fastq ../$ref.$j $readPair ../$lib1 $lengthlib1 $insertlib1 ../$lib2 $lengthlib2 $insertlib2 ../$ref.$i.shift.txt 3 $i $minQual"

		bsub -q normal  -w"$nameSecondStage" -J "Clean.$tmpId.$root$i" -e Clean.e -o Clena.o "$ICORN_HOME/clean.ssahapileup.sh $root $num"
	else
	    $ICORN_HOME/pileup.joinData.sh $root ../$root.mate.fastq ../$ref.$i 0 $genomesize $readPair $insertSize &> tmp.pileup.$i.out
	    $ICORN_HOME/helper2.sh $root
echo " perl $ICORN_HOME/icorn.Correct.pl ../$ref.$i $root.snp  $root.del $root.ins Ins.Reads.fastq ../$ref.$j $readPair ../$lib1 $lengthlib1 $insertlib1 ../$lib2 $lengthlib1 $insertlib1 ../$ref.$i.shift.txt 3 $i $minQual"
	    perl $ICORN_HOME/icorn.Correct.pl ../$ref.$i $root.snp  $root.del $root.ins Ins.Reads.fastq ../$ref.$j $readPair ../$lib1 $lengthlib1 $insertlib1 ../$lib2 $lengthlib1 $insertlib1 ../$ref.$i.shift.txt 3 $i $minQual

	

		$ICORN_HOME/clean.ssahapileup.sh $root 1
	fi #end if LSF=1
	
	cd ..

	## get the result lines
	$ICORN_HOME/icorn.CollectResults.sh $ref $root
	export RefAll; grep '>' $RefAll | sed 's/>//g' | sed 's/\s+//g' | perl -nle '($contig) = $_ =~ /(\S+)/; system "cat $ENV{RefAll}.*.$contig.gff > All.$contig.gff"'
done;