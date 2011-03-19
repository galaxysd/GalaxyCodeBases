#!/usr/bin/perl-w
use strict;
#explanation:this program is edited to 
#edit by HeWeiMing;
 die  "Version 1.0 2009-7-6;\nUsage: $0 <InPut_RST><InPut_addc><OutDir>\n" unless (@ARGV == 3);

open	 A,"$ARGV[0]"  || die "$!" ;
open     B,"$ARGV[1]"  || die "$!" ;
open     C,">$ARGV[2]"  || die "$!" ;


my %hash=();
	while(<A>) 
	{ 
		chomp ; 
		my @inf=split ;
		$hash{$inf[0]."_".$inf[1]}=$inf[2]	;
	}
close A;

	while(<B>)
	{
		chomp ;
        	my @inf=split ;
		my $line1=join("\t",@inf[0..$#inf]);
		print  C  $line1,"\t",$hash{$inf[0]."_".$inf[1]},"\n";
	}

close B;
close C;




#		my %iupac = ("M", "ac", "K", "gt", "Y", "ct", "R", "ag", "W", "at", "S", "cg","A","aa","T","tt","C","cc","G","gg","-","--");
#		foreach my $ke (keys %iupac)
#		{
#		$inf=~s/$ke/$iupac{$ke}/g;
#		}
#						$a1[0]=0+($inf=~s/a/A/g);	 
#						$a1[1]=0+($inf=~s/t/T/g);	
#						$a1[2]=0+($inf=~s/c/C/g);	 
#						$a1[3]=0+($inf=~s/g/G/g);	
#                    ($aSed)=(sort{$a<=>$b} @a1)[2];
#foreach my $k (sort {($a <=> $b )or( $hash{$a} <=> $hash{$b}) }keys %hash)
#	{
#	print B	$k,"\t", $hash{$k},"\n";
#	}
#close B;
#
#		chomp ; 
#		my @inf=split/\t/ ;
#		my $line=join("\t", @inf[0..$#inf]);

#my $file=`ls /share/*/*.maf`;
#my @arry=split/\n/, $file;
#
#
#my $leng_file=@arry;
#for (my $ii=0;$ii<$leng_file;$ii++)
#{
#                print $arry[$ii],"\n";
#
#        open     A,"$arry[$ii]"  || die "$!" ;
#        
#}
