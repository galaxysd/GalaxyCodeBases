#!/usr/bin/env perl
#############################################
#Author: Jiang Li
#email: riverlee2008@gmail.com
#Creat Time: Fri 03 Aug 2012 10:34:54 PM CDT  
#Vanderbilt Center for Quantitative Sciences
#############################################
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::File;

######### Update from pile2base.pl ####################
# 1).it will parse the insert and deletion too, 
# 2).it can parse pileups from multiple sample,
# 3).you don't need to define the outputfile,instead, 
#    you will need to provide 'prefix' which will name result from each sample
#    as prefix1.txt, prefix2.txt,etc. Default, prefix is sample.
# 4).Provide option to define the base quality score offset, default is 33 (Sanger standard).
# 5).Read parameters from command line with options

#Usage: perl pileup2baseindel.pl -i <pileupfile> -bq [BQcutoff] -prefix [sample] -offset [33]
my $usage = <<USAGE;
Usage: perl pileup2base.pl -i <pileupfile> -bq [BQcutoff] -prefix [sample] -offset [33]
        -i        input pileup file, could be from 1 sample or multiple samples
        -s        pileupfile with mapping quality
        -bq       base quality score cutoff for each mapped/unmapped base,
                  only those larger than cutoff will be output in the result, default is -5, means no filter
        -prefix   output file prefix, default is sample, the output will be named as prefix1.txt, prefix2.txt, etc.
        -offset   Offset to change ASCII character to base quality score, default is 33 (sanger format).
        -d        DUMP mode
        -h        print out this
USAGE

my ($input,$BQcut,$offset,$prefix,$help,$hasMQ,$doDUMP) = (undef,-5,33,"sample",undef,undef);
GetOptions(
	"i=s"=>\$input,
	"bq=i"=>\$BQcut,
	"offset=i"=>\$offset,
	"prefix=s"=>\$prefix,
	"s"=>\$hasMQ,
	"d"=>\$doDUMP,
	"h"=>\$help
);

if($help){
	print $usage;
	exit(0);
}
my $itemLen = 3;
$itemLen = 4 if $hasMQ;

unless ($input){
	print "Input file does not provide yet\n";
    print "\n$usage\n";
    exit(1);
}
if(! -e $input){
    print "Input file '$input' does not exists\n";
    print "\n$usage\n";
    exit(1);
}

#Do the parsing
open FILE, $input or die "error, can not open $input";
print "[",scalar(localtime),"] Begin parsing...\n";

my $line = <FILE>;
$line=~s/\r|\n//g;
my ($chr,$loc,$ref,@dp_bases_bq) = split /\s+/, $line;
my $n = int (scalar(@dp_bases_bq)/$itemLen); #determine how many samples, use int just in safe
my %files;
foreach my $i (1..$n){
	my $fh = new IO::File;
	$fh->open("> ${prefix}${i}.txt");
	if ($doDUMP) {
		print $fh join("\t",qw(chr loc ref Depth La Lc Lg Lt Bases BaseQs Insertion Deletion)),"\n";
	} else {
		print $fh "chr\t"."loc\t"."ref\t"."A\t"."T\t"."C\t"."G\t"."a\t"."t\t"."c\t"."g\t"."Insertion\t"."Deletion\n";
	}
	$files{$i}=$fh;
}

my $theFunc = \&parsePileup;
$theFunc = \&dumpPileup if $doDUMP;

seek FILE, 0, 0;
while(<FILE>){
	s/\r|\n//g;
	my ($chr,$loc,$ref,@dp_bases_bq) = split /\s+/;
	my $n = int (scalar(@dp_bases_bq)/$itemLen); #determine how many samples, use int just in safe
	foreach my $i (1..$n){
		my $fh = $files{$i};
		my @region=($itemLen*($i-1),$itemLen*($i-1)+1,$itemLen*($i-1)+2);
		my ($dp,$bases,$bq) = @dp_bases_bq[@region];
		my $str = $theFunc->($ref,$bases,$bq,$BQcut,$offset,$dp);
		if($str ne "*"){
			print $fh join "\t",($chr,$loc,$ref,$str);
		}
	}
}

close FILE;
foreach my $k (keys %files){
	$files{$k}->close;
}
print "[",scalar(localtime),"] Finished\n";

sub dumpPileup($$$$$$) {
	my ($ref,$bases,$bq,$BQcut,$offset,$dp) = @_;
	if($bases eq "*"){
		return "*";
	}
	$bases=~s/\^.//g;
	$bases=~s/\$//g;
	my %hash=();
	my %deletion=();
	while($bases=~/-(\d+)/g){
		$hash{$1}=1;
	}
	foreach my $k (keys %hash){
		while($bases=~/-$k([ACGTNacgtn]{$k})/g){
			$deletion{$1}++;
		}
		$bases=~s/-$k[ACGTNacgtn]{$k}//g;
	}
	
	%hash=();
	my %insertion=();
	while($bases=~/\+(\d+)/g){
		$hash{$1}=1;
	}
	foreach my $k (keys %hash){
		while($bases=~/\+$k([ACGTNacgtn]{$k})/g){
			$insertion{$1}++;
		}
		$bases=~s/\+$k[ACGTNacgtn]{$k}//g;
	}
	my @base=split (//,$bases);
	my @bq=split(//,$bq);
	my (@nq,@newBases);
	#my $realdep=0;
	my $skipped=0;
	for(my $i=0;$i<@base;$i++){
		my $ch=$base[$i];
		if($ch eq "."){
			push @newBases,uc $ref;
		}elsif($ch eq ","){
			push @newBases,lc $ref;
		}elsif($ch eq '*' or $ch eq '#'){
			++$skipped;
			next;
		}else{
			push @newBases,$ch;
		}
		my $score=ord($bq[$i])-$offset;
		push @nq,$score;
		#++$realdep;
	}
	#if(scalar(@newBases) ne scalar(@nq)){
	#	die ">>> $_\n";
	#}
	my $sbases = join('',@newBases);
	my $sbq = join(',',@nq);

	my $insertion="NA";
	my $deletion="NA";
	if(scalar(keys %insertion)){
		$insertion="";
		foreach my $k (sort {$insertion{$b}<=>$insertion{$a} || $b cmp $a} keys %insertion){
			$insertion.=$insertion{$k}.":".$k."|";
		}
		chop($insertion);
	}
	
	if(scalar(keys %deletion)){
		$deletion="";
		foreach my $k (sort {$deletion{$b}<=>$deletion{$a} || $b cmp $a} keys %deletion){
			$deletion.=$deletion{$k}.":".$k."|";
		}
		chop($deletion);
	}
	my ($La,$Lc,$Lt,$Lg)=(0,0,0,0);
	for(my $i=0;$i<@newBases;$i++){
		my $ch=$newBases[$i];
		my $score=$nq[$i];
		if ($score>40) {
			#print STDERR '^';
			$score=40;
		}
		my $Erate = 10**(-$score/10);
		#print "$ch $score $Erate\n";
		if ($ch eq 'A' or $ch eq 'a') {
			$La += log(3*(1-$Erate)/$Erate);
		} elsif ($ch eq 'C' or $ch eq 'c') {
			$Lc += log(3*(1-$Erate)/$Erate);
		} elsif ($ch eq 'T' or $ch eq 't') {
			$Lt += log(3*(1-$Erate)/$Erate);
		} elsif ($ch eq 'G' or $ch eq 'g') {
			$Lg += log(3*(1-$Erate)/$Erate);
		} else {
			die "[$ch] @newBases\n$ref,$bases,$bq,$BQcut,$offset,$dp\n";
		}
		#print "$La,$Lc,$Lt,$Lg\n";
	}
	my $str=join("\t",$dp-$skipped,$La,$Lc,$Lg,$Lt,$sbases,$sbq,$insertion,$deletion)."\n";
	return $str;
}

sub parsePileup{
	my ($ref,$bases,$bq,$BQcut,$offset) = @_;
	
	if($bases eq "*"){
		return "*";
	}
	
	#do some modificaton on $base to remove additional characters
	#1,remove the ^. pattern
	$bases=~s/\^.//g;
	#2,remove the $ pattern
	$bases=~s/\$//g;
	#3,remove -[0-9]+[ACGTNacgtn]+ pattern
	my %hash=();
	my %deletion=();
	while($bases=~/-(\d+)/g){
		$hash{$1}=1;
	}
	#get the deletion sequences and delete them
	foreach my $k (keys %hash){
		while($bases=~/-$k([ACGTNacgtn]{$k})/g){
			$deletion{$1}++;
		}
		$bases=~s/-$k[ACGTNacgtn]{$k}//g;
	}
	
	%hash=();
	my %insertion=();
	while($bases=~/\+(\d+)/g){
		$hash{$1}=1;
	}
	foreach my $k (keys %hash){
		while($bases=~/\+$k([ACGTNacgtn]{$k})/g){
			$insertion{$1}++;
		}
		$bases=~s/\+$k[ACGTNacgtn]{$k}//g;
	}

	#Now @base and @bq have the same length
	my @base=split (//,$bases);
	my @bq=split(//,$bq);
	#I have check it
	#if(scalar(@base) ne scalar(@bq)){
	#	print $_,"\n";
	#}
	#foreach my $c (@base){
	#	$check{$c}++;
	#}
	my $forward_A=0;
	my $forward_T=0;
	my $forward_C=0;
	my $forward_G=0;
	my $reverse_A=0;
	my $reverse_T=0;
	my $reverse_C=0;
	my $reverse_G=0;

	#start the loop
	for(my $i=0;$i<@base;$i++){
		my $ch=$base[$i];
		my $score=ord($bq[$i])-$offset;
		if($score>=$BQcut){
		    if($ch eq "A"){
		        $forward_A++;
		    }elsif($ch eq "T"){
		        $forward_T++;
		    }elsif($ch eq "C"){
		        $forward_C++;
		    }elsif($ch eq "G"){
		        $forward_G++;
		    }elsif($ch eq "a"){
		        $reverse_A++;
		    }elsif($ch eq "t"){
		        $reverse_T++;
		    }elsif($ch eq "c"){
		        $reverse_C++;
		    }elsif($ch eq "g"){
		        $reverse_G++;
		    }elsif($ch eq "."){
		        if($ref eq "A"){
		            $forward_A++;
		        }elsif($ref eq "T"){
		            $forward_T++;
		        }elsif($ref eq "C"){
		            $forward_C++;
		        }elsif($ref eq "G"){
		            $forward_G++;
		        }
		    }elsif($ch eq ","){
		        if($ref eq "A"){
		            $reverse_A++;
		        }elsif($ref eq "T"){
		            $reverse_T++;
		        }elsif($ref eq "C"){
		            $reverse_C++;
		        }elsif($ref eq "G"){
		            $reverse_G++;
		        }
		    }
		}#end the condition  $score>=$BQcut
	}#end the loop
	
    #my $str="$chr\t$loc"."\t".$ref."\t".$forward_A."\t".$forward_T."\t".$forward_C."\t".$forward_G."\t".$reverse_A."\t".$reverse_T."\t".$reverse_C."\t".$reverse_G."\t";
	my $str=$forward_A."\t".$forward_T."\t".$forward_C."\t".$forward_G."\t".$reverse_A."\t".$reverse_T."\t".$reverse_C."\t".$reverse_G."\t";
	my $insertion="NA";
	my $deletion="NA";
	if(scalar(keys %insertion)){
		$insertion="";
		foreach my $k (sort {$insertion{$b}<=>$insertion{$a} || $b cmp $a} keys %insertion){
			$insertion.=$insertion{$k}.":".$k."|";
		}
		chop($insertion);
	}
	
	if(scalar(keys %deletion)){
		$deletion="";
		foreach my $k (sort {$deletion{$b}<=>$deletion{$a} || $b cmp $a} keys %deletion){
			$deletion.=$deletion{$k}.":".$k."|";
		}
		chop($deletion);
	}
	$str.=$insertion."\t".$deletion."\n";
	return $str;
}

