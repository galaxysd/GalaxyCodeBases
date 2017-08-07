#!/bin/bash

res=$1;
fastq=$2;
ref=$3;
trans=$4;
genome_size=$(seqstat $ref | awk '/Total/ {print $4}')
readPair=$5;
InsertSize=$6
givenCigar=$7

function printit {
##	echo `echo $1 | awk -v FS=" " '{print $1}'` >> $res.txt
	echo $1 >> $res.txt
}

echo " The ssaha pileup script should have been used. IMPORTANT, to just have data of one run in this directory!\nStatistical result is written in $res.txt";
echo $givenCigar


### get the insert size
minInsertSize=$(echo $InsertSize | perl -nle '/^(\d+),(\d+)$/;print $1')
maxInsertSize=$(echo $InsertSize | perl -nle '/^(\d+),(\d+)$/;print $2')

if [ -z "$givenCigar" ]; then
	if [ "$readPair" = "1" ]; then
        cat *cigar3* > $res.cigarX
		
		perl  $ICORN_HOME/analyse.matePairDistance.RP.pl $res.cigarX $res $maxInsertSize $minInsertSize
	else
		cat *cigar3* >  $res.cigar
	fi;
	cat *cigar1* > $res.AllMapped.cigar
else
	res=$givenCigar;
	echo "hello";
fi;

prog=$PILEUP_HOME/ssaha_pileup;





### name to file
echo "$res" > $res.txt
echo " " >> $res.txt

### genome size
echo $genome_size >> $res.txt
echo " " >> $res.txt

### get amount of reads in fastq file $fastq
out=$(wc $fastq | awk '{print $1}')
reads=$(($out/4));
echo "$reads" >> $res.txt


# amount of read that map at all
echo "Amount of reads map in general, also multiple";
out=$(awk '{ print $2 }'  $res.AllMapped.cigar | uniq | wc | awk '{print $1}')
echo  $out >> $res.txt;

num=$(($out*100/$reads));
echo "$num" >> $res.txt; 
echo " " >> $res.txt


### amount of read that map uniquely
echo "Amout of reads, after filtering";
out=$(awk '{ print $2 }' $res.cigar | uniq | wc| awk '{print $1}')
printit $out;
num=$(($out*100/$reads));
echo "$num" >> $res.txt
#printit $num; 
echo " " >> $res.txt

if [ "$readPair" = "1" ]; then
#### echo "amount pairs";
	perl $ICORN_HOME/analyse.countpair.pl  $res.cigar $maxInsertSize $minInsertSize  >> $res.txt
else
	echo >> $res.txt;
	echo >> $res.txt;
	echo >> $res.txt;
	echo >> $res.txt;
fi

# get the correct fastq order!
awk '{print $2}' $res.cigar > $res.fastq.names
$prog/other_codes/get_seqreads/get_seqreads $res.fastq.names $fastq $res.fastq


$prog/ssaha_pileup/ssaha_pileup -solexa 1 -trans $trans $res.cigar $ref $res.fastq > $res.snp

$prog/ssaha_pileup/ssaha_pileup -cons 1 -trans $trans -solexa 1 $res.cigar $ref $res.fastq > $res.pileup
$prog/ssaha_pileup/ssaha_indel -insertion 1 $res.cigar $ref  $res.pileup > $res.ins
$prog/ssaha_pileup/ssaha_indel -deletion 1  $res.cigar $ref  $res.pileup > $res.del


if [ "$SsahaDoCoveragePlots" = "1" ] ; then
	if [ ! -d plot ]		# be sure the directory /mnt exists
		then
		mkdir plot
		
	fi
	awk '{ if ($1 ~ "^cons") {print $4 > "plot/"x"."$2".plot"}}' x=$res $res.pileup
fi

#perl ~tdo/Bin/pileup.parsePileup2Plot.pl $res.pileup $res
echo " " >> $res.txt


### echo "None covered bp of genome";
out=$(awk '{if ($4==0){print $4} }' $res.pileup | wc | awk '{print $1}')
cov=$(($out*100/$genome_size))
printit $out;
echo $cov >> $res.txt

echo "Covered < 5 bp of genome";
out=$(awk '{if ($4<5){print $4} }' $res.pileup | wc | awk '{print $1}')
cov=$(($out*100/$genome_size))
printit $out;
echo $cov >> $res.txt

### echo "Covered < 10 bp of genome";
out=$(awk '{if ($4 < 10){print $4} }' $res.pileup | wc | awk '{print $1}')
cov=$(($out*100/$genome_size))
printit $out;
echo $cov >> $res.txt

echo "Covered >= 10 bp of genome";
out=$(awk '{if ($4 >= 10){print $4} }' $res.pileup | wc | awk '{print $1}')
cov=$(($out*100/$genome_size))
printit $out;
echo $cov >> $res.txt

echo "Covered >= 20 bp of genome";
out=$(awk '{if ($4 >= 20){print $4} }' $res.pileup | wc | awk '{print $1}')
cov=$(($out*100/$genome_size))
printit $out;
echo $cov >> $res.txt

echo "Covered >= 40 bp of genome";
out=$(awk '{if ($4 >= 40){print $4} }' $res.pileup | wc | awk '{print $1}')
cov=$(($out*100/$genome_size))
printit $out;
echo $cov >> $res.txt


echo " " >> $res.txt
echo " " >> $res.txt
echo "Amount SNPs Qual >= 20";
out=$(awk '{if ($3 >= 20){print $4} }' $res.snp | wc)
printit $out;

echo "Amount SNPs Qual >= 50";
out=$(awk '{if ($3 >= 50){print $4} }' $res.snp | wc)
printit $out;

echo "Amount SNPs Qual >= 60";
out=$(awk '{if ($3 >= 60){print $4} }' $res.snp | wc)
printit $out;

echo "Amount SNPs Qual >= 80";
out=$(awk '{if ($3 >= 80){print $4} }' $res.snp | wc)
printit $out;

echo "Amount SNPs Qual >= 90";
out=$(awk '{if ($3 >= 90){print $4} }' $res.snp | wc)
printit $out;

