#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use Getopt::Long;
use FindBin qw($Bin);

#Edit by heweiming
#Bulit 2009-08-01
#last modified 2009-08-14

#this program calculate the copy_number,5bp_edge_reads,No_5bp_edge_reads,Uniq_reads,ration_check,RSTcheck,and filter the result by the rule as :
#1.depth [20-250], 2.quality>=15 ,3.copy_number<=1.5
#4.there is at least one read that not 5bp_edge coveraged this site.
#5.at least one reads uniq corveraged this site
#6.ration_check<=2
#7.RST_check>=0.005

##note:the rule is only for rice,if you do other case,you must change the filtering rule .

## format_of_copyNumRST5p_v2.pl.output
#chr loc refBase depth cpnum base1 base2 numOfBase1 numOfBase2 quality copynum 5bp_edge no_edge Uniq_reads ration_check RST_check
#chromosome01    173     C       95      0.594737        C       T       98      4       18      1.17    27      68      84      0.00    0.364752


die  "Version 1.0 2009-7-6;\nUsage:
perl $0 -i population snp  -r reference Chr01.fa -l chr_length -c Chr -n 50 -m soap.list -o output\n" unless (@ARGV == 12);
# copyNumLcRst.pl -i ./Watermelon_17/SNP/PE/seg01/seg01 -r ./watermelon_v2_888/faByChr/seg01.fa -l ./Watermelon_17/FinalSNP2/watermelon.merge.len -c seg01 -n 16 -m ./Watermelon_17/SortbyChr/PE/List/seg01.list -o ./Watermelon_17/FinalSNP2/FinalSNP/Population/seg01

my ($numberOfFile,$reference,$length_chr_file,$chromosome,$input,$mergelist,$outfile,$help);
GetOptions(
#	"numberOfFile:i"=>\$numberOfFile,	# Sample count
	"chromosome:s"=>\$chromosome,	# ChrID
	"length_chr_file:s"=>\$length_chr_file,	# ChrLen or Chr.nfo
	"reference:s"=>\$reference,
	"input:s"=>\$input,	# pCNS
	"merge:s"=>\$mergelist,	# merged soap list, by Chr. Changed to accept univsrtal megred.lst
	"out:s"=>\$outfile,
	"help"=>\$help,
);
#$numberOfFile ||= 1;
$chromosome ||= "Gm01";

open(RAW,"$input") || die"$!";  # open my raw soapsnp data
my %exis_snp=();

while(<RAW>)
{
	chomp;
        my @line=split(/\s+/,$_);
	my $inf=join("\t", @line[0..$#line]);
        if($line[3]>= 1  && $line[9]>= 1 ) 		# depth > 0   and quality > 0
	{  $exis_snp{$line[0]."_".$line[1]}=$inf;}
}
close RAW;

my %hash;
my %noGapLength;
open L,$length_chr_file || die $!;
while (<L>){
	chomp;
	next if (/total/);
	my @words = split;
	$hash{$words[0]} = $words[1];
#	$noGapLength{$words[0]} = $words[1] - $words[2];
}
close L;

$_="";


open A,'<',$mergelist || die "$!" ;
my @cpnumber;
my @hit;
my %ina;
my %outa;
my %hashLC=();
my %hashRST=();
my %uniqHash=();
my $length = $hash{$chromosome};



# mkdir $depthpath unless -e $depthpath ;

while(my $line=<A>)
{
	chomp $line;
	my ($Sample,$Chr,$Len,$file)=split /\t/,$line;
	next if $Chr ne $chromosome;

	open IN,'<',$file or die $!;

	while (<IN>){
		chomp;
		my ($word1,$word2,$hit,$read_len,$chr,$pos)=(split /\s+/)[1,2,3,5,7,8];
		next unless $chr eq $chromosome;
		my $start = $pos;
		my $end = $pos + $read_len - 1;
	#	if($start>$print_end+1)
	#	{
	#		my $tempddd="";
	#		for(my $jj=$print_end; $jj<$start; $jj++)
	#		{
	#		if (!exists $locadethp[$jj]) {$locadethp[$jj]=0;}
	#		#print  OUTSoapdethp  $locadethp[$jj]," ";
	#		$tempddd.=$locadethp[$jj]." ";
	#		$locadethp[$jj]=();  delete   $locadethp[$jj] ;
	#		}
	#		$print_end=$start;
	#		print OUTSoapdethp  $tempddd;
	#		$tempddd=();
	#	}

		foreach my $j ($start..$end){
			my $loc_chr_a=$chr."_".$j;
	#		$locadethp[$j]++;
			if(exists $exis_snp{$loc_chr_a})
			{
				$hit[$j] += $hit;
				$cpnumber[$j] ++;

				if($j-$start<=5  ||  $end-$j<=5)
				{
				$ina{$loc_chr_a}++;
				}
				else
				{
				$outa{$loc_chr_a}++;
				}

			 	next if ($hit != 1) ;
				$uniqHash{$loc_chr_a}++;
				my $offset = $j - $start;
				my $allele = substr($word1,$offset,1);
				my $q = substr($word2,$offset,1);
				$hashRST{$loc_chr_a}.= $allele.$q;
				$hashLC{$loc_chr_a}.= $allele.$q;
			}
		}
	}
	close IN;
	foreach (keys %exis_snp){
		$hashLC{$_} .= "~,";
	}



#	close  OUTSoapdethp ;
#	 $outputfile=();
#	  @locadethp=();
 #        $print_end=();
}

close A;
$_="";
# $depthpath=();
 print "read soap had done!\n";
##sereach position of gap
my %ref;
open (A,$reference) || die $!;
$/=">";
<A>;
while(<A>){
    chomp ;
    my @inff=split/\n/;
	my $first=shift @inff;
	my $id= (split/\s+/,$first)[0];
	$ref{$id}=join("",@inff[0..$#inff]);
	@inff=();

#	chomp;
#	if(/^>(\S+)/){
#		$id = $1;
#	}else{
#		$ref{$id} .= $_;
#	}
}

$/="\n";

$_="";

close A;
$/="\n";


foreach my $key (keys %ref){
	next if ($key ne $chromosome);
	my $index = length ($ref{$key});
	foreach my $i (0..$index-1){
		if (substr($ref{$key},$i,1) eq "N"){
			my $position = $i + 1;
            if(exists $cpnumber[$position])
            {
			$cpnumber[$position] = -1;
            }
		}
	}
	$ref{$key}="";
}
##sereach position of gap
%ref=();
undef %ref;

$_="";
# delete %ref;


 # my $sum = 0;
 # my $totalnoGapLength = 0;

foreach my $h(1..$length){
if(	!exists $cpnumber[$h] )
{
    next ;
}
	if ($cpnumber[$h] != 0 && $cpnumber[$h] != -1){
		$cpnumber[$h] = sprintf ("%.2f",$hit[$h] / $cpnumber[$h]);
#		$sum++;
	}
#	if ($cpnumber[$h] != -1)
#	{
#		$totalnoGapLength++;
#	}
	$hit[$h]=();
	delete $hit[$h];
}

@hit=();
undef @hit;
$_="";
# delete @hit;


#my $coverage = $sum / $totalnoGapLength;

my  $outDir=$outfile;

#   open(Report,">$outDir\.report")||die"$!";

#print Report  "chr_name\tchr_len\tchr_noGapLength\ttotalnoGapLength\tsum\tcoverage\n";
#print Report  "$chromosome","\t",$hash{$chromosome},"\t",$noGapLength{$chromosome},"\t$totalnoGapLength\t$sum\t$coverage\n";

#close  Report ;

#----- got each copynumber  and filter
#my @cp =  @t[1..$length];

# @t=();
  %noGapLength=();  %hash=(); $length=();
   undef  %noGapLength ; undef  %hash ; undef $length ;

 #delete $sum;
     open(Addcn,">${outDir}.add")||die"$!";  # print add_cn out
     open(RST,">${outDir}.rs")||die"$!"; # print filter SNP
 #    open(OUTdb,">$outDir\.dbsnp")||die"$!"; # print filter SNP


my %allele = ();
my %count = ();


#foreach my $k (sort   keys %exis_snp)
foreach my $k (sort {my ($a,$b),(split(/\_/,$a))[1]<=>(split(/\_/,$b))[1]}keys %exis_snp)
	{
	my $inf=$exis_snp{$k};
	my @line=split(/\s+/,$inf);
	my $copynum = $cpnumber[$line[1]];
	my $LC_input=$hashLC{$k};
	my $lineprint=();
	 $hashLC{$k}=();
 $cpnumber[$line[1]]=();
# delete   $cp[$line[1]-1];

#	chomp $inf;           ADDcn  OUT
if(!exists  $ina{$k} ) {$ina{$k}=0;}
if(!exists $outa{$k} ) { $outa{$k}=0;}
if(!exists  $uniqHash{$k}) { $uniqHash{$k}=0; }
#  print Addcn  $inf,"\t",$copynum,"\t",$ina{$k},"\t",$outa{$k},"\t", $uniqHash{$k};
 $lineprint= $inf."\t".$copynum."\t".$ina{$k}."\t".$outa{$k}."\t". $uniqHash{$k};
delete  $ina{$k} ; delete  $outa{$k} ;  $exis_snp{$k}=0; delete  $uniqHash{$k};
delete $hashLC{$k};
 undef   $ina{$k} ;  undef   $outa{$k}  ;  undef   $uniqHash{$k} ; undef     $hashLC{$k};
#			LC   ratio  OUT
  my @tt= split (/\,/,$LC_input);  ##if not covered, it will report error
  $LC_input=();
	for (my $j=0;$j<=$#tt;$j++)
	{
                my @a = split(//,$tt[$j]);
                my $seq = "";
                my $quality = "";
                for (my $k=0;$k<=$#a;$k+=2){
#                       next if ($a[$k] eq "~");
                        $seq .= $a[$k];
                        if ($#a>0 && $a[$k] ne "~"){
                                $quality .= $a[$k+1];
                        }
                }
                $LC_input .= $seq.$quality.",";
        }


 my @snp = (split(/\t/,$inf))[5,6];
        my %snp;
        foreach (@snp){
                $snp{$_} = 1;
        }
        my %stat;my %het=();
        my @reads = split(/\,/, $LC_input);
#       print $key,"\t"
        for (my $i=0;$i<@reads;$i++){
                if ($reads[$i] =~ /^~$/){
#                       print "- ";
                }else{
                    my $seq = (split(/\~/,$reads[$i]))[0];my $quality = (split(/\~/,$reads[$i]))[1];
                        my @allele_LC= split(//,$seq);my @q=split(//,$quality);
                        for (my $k = 0;$k<@q;$k++)
                        {
                         $q[$k] = ord($q[$k]) - 64;
                        }
                        my %hashtemp;
                        for (my $j = 0;$j<@allele_LC;$j++){
                                next if (! exists $snp{$allele_LC[$j]});
                                next if ($q[$j]<5);
                                $hashtemp{$allele_LC[$j]} ++;
                                $stat{$allele_LC[$j]} ++;
                        }
                        my $l = scalar keys %hashtemp;

			  if ($l == 1){
                                foreach (keys %hashtemp){
#                                       print $hash{$_},$_;
                                }
                        }
                        if ($l > 1){
                                my $output = "";
                                foreach (sort keys %hashtemp){
                                        $het{$_} += $hashtemp{$_} if ($l > 1);
                                        $output .= $hashtemp{$_}.$_.":";
#                                       print $hash{$_},$_,":";
                                }
                                chop $output;
#                               print $output;
                        }
#                       print " ";
                }
        }
        my @just;my @stat_reverse;
        foreach (sort keys %stat){
##              print $stat{$_},$_,":";
                push @just,[$stat{$_},$_];
        }
        @just = sort{my ($a,$b),$b->[0]<=>$a->[0] || $a->[1] cmp $b->[1]} @just;
#        print OUT $just[0][0],$just[0][1],":",$just[1][0],$just[1][1];
#       print OUT "\t";
       if (scalar keys %het > 1){
               # print OUT  $het{$just[0][1]},$just[0][1],":",$just[1][0],$just[1][1];
                my $rate = sprintf ("%.2f", $het{$just[0][1]}/$just[1][0]);
               #  print Addcn "\t$rate";
		$lineprint=$lineprint."\t".$rate."\n";
#               my $hom = $just[1][0] - $het{$just[1][1]};
#               my $hom = $just[1][0] - $het{$just[1][1]};
#               my $rate2 = $hom / $het{$just[1][1]};
#               print "\t$rate2";
        }else{
          #      print OUT "0:0\t0";
#                print Addcn "\t",0;
		 $lineprint=$lineprint."\t0.00\n";
        }

        print Addcn   $lineprint;
		$lineprint=();  undef    $lineprint;



@tt=();
#		RST  out
	  @tt = split (//,$hashRST{$k}); ##if not covered, it will report error
                for (my $jl=0;$jl<$#tt;$jl++){
                        if ($jl % 2 == 0){
                                $allele{$tt[$jl]} .= $tt[$jl+1];
                                $count{$tt[$jl]} ++;
                        }
                }

                my $temp;
                if (exists $allele{'A'}){
                        $temp = $allele{'A'} . "\t";
                }else{
                        $temp = "." . "\t";
                }
                if (exists $allele{'C'}){
                        $temp .= $allele{'C'} . "\t";
                }else{
                        $temp .= "." . "\t";
                }
                if (exists $allele{'T'}){
                        $temp .= $allele{'T'} . "\t";
                }else{
                        $temp .= "." . "\t";
                }
                if (exists $allele{'G'}){
                        $temp .= $allele{'G'};
                }else{
                        $temp .= ".";
                }
                %allele = ();
                %count = ();
                $hashRST{$k} = $temp;
                print RST  $inf,"\t",$copynum,"\t",$hashRST{$k},"\n";
		 $hashRST{$k}=();
		delete   $hashRST{$k} ;
		undef    $hashRST{$k} ;
	}

close Addcn;
close RST;

%hashRST=();   undef %hashRST;
%ina=();	undef %ina;
%outa=();	undef %outa;
%hashLC=();	undef %hashLC;


#delete %ina;
#delete  %outa;
#delete %hashRST;

my $file1="${outDir}.add";
my $file2="${outDir}.rs";
my $file4="${outDir}.add_cn";
my $file3="${outDir}.rst";
system ("$Bin/RankSumTest  -i   $file2  -o   $file3 ");
system ("perl $Bin/addRst.pl  $file3   $file1  $file4 ");
#unlink    $file2 ;          #   system("rm  $file2");
#unlink    $file1 ;          #   system("rm $file1");
#unlink    $file3 ;          #   system("rm $file3");

###############  filter the addcn ,you must modif it since the level of filter ###########
#system("perl  /panfs/GAG/junli/raid010/pipeline/populationsnp/bin/filter_addcn.pl  $file4 "  );

