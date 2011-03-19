#!/bin/env perl
use strict;
use warnings;
# original from HeWeiMing;

die "Version 1.1 2010-7-22;\nUsage:$0 <sample number> <indivi_snp file> <output prefix> <Isselfedseed> <cutoff>\n" unless (@ARGV >= 1 && @ARGV < 5);

open B,'<',$ARGV[1]  || die "$!" ;
open C,'>',"$ARGV[2]\.add"  || die "$!" ;
open D,'>',"$ARGV[2]\.add.filter"  || die "$!" ;
open E,'>',"$ARGV[2]\.rm"  || die "$!" ;
my $Isselfedseed=1 ;
my $cuoff=11.34487 ;
$ARGV[0]||=33;
my $SampleCount=$ARGV[0];
if (exists $ARGV[3])
{
    $Isselfedseed= 0 ;
}
if(exists $ARGV[4])
{
    if ( $ARGV[4] <1 &&  $ARGV[4] >0.995) {$cuoff=10.59663 ;}
    elsif( $ARGV[4] >=0.975 &&  $ARGV[4] <0.99) {$cuoff=7.37776 ;}
    elsif( $ARGV[4] >=0.95 &&  $ARGV[4] <0.975) {$cuoff= 5.99146 ;}
    elsif( $ARGV[4] >=0.90 &&  $ARGV[4] <0.95 ) {$cuoff= 4.60517 ;}
}

my $header=<B>;
if ($header =~ /^#/) {
	my @Samples=split /,\s*/,(split /\t/,$header)[1];
	$header="#ChrID\tPos\t".join(' ',@Samples);
} else {
	seek B,0,0;
	$header=undef;
}
if (defined $header) {
	print C $header,"\tValues\n";
	print D $header,"\n";
}
while(<B>)
{
	chomp;
	my @snp=split /\s+/;
	my @inf=@snp[0,1];
	my @gtypes;
	for(my $i=3;$i<=$#snp-3;$i+=3) {
		if ($snp[$i+2]==0) {$snp[$i]="-";}
		push @gtypes,$snp[$i];
	}
	push @inf,join(' ',@gtypes);
	#my @inf=split /\t/;
	my @a1 ;
	my $yy=$inf[2];
	$a1[0]=0+($yy=~s/A/1/g);
	$a1[1]=0+($yy=~s/C/2/g);
	$a1[2]=0+($yy=~s/T/3/g);
	$a1[3]=0+($yy=~s/G/4/g);
	my $e=0+($yy=~s/-/5/g);
#B observed A theor;
	my $BB=(sort{$a<=>$b} @a1)[3];
#	print $BB,"\t";
	my $bb=(sort{$a<=>$b} @a1)[2];
#	print $bb,"\t";
	my $Bb=$SampleCount-$e-$BB-$bb;
	my $sum=$SampleCount-$e;
	my $p=(2*$BB+$Bb)/(2*$sum);
	my $q=1-$p;
#	my $sum=($BB+$Bb+$bb);
	my $AA=$p*$p*$sum;
	my $aa=$q*$q*$sum;
	my $Aa=2*$q*$p*$sum;
	my $p1; my $p2; my $p3;
	if($AA!=0){ $p1=((abs($BB-$AA))*(abs($BB-$AA))/$AA);} else {$p1=0;}
	if($Aa!=0){ $p2=((abs($Bb-$Aa))*(abs($Bb-$Aa))/$Aa);} else {$p2=0;}
	if($aa!=0){ $p3=((abs($bb-$aa))*(abs($bb-$aa))/$aa); } else {$p3=0;}
	my $x2=$p1+$p2+$p3;
#	print $AA,"\t",$Aa,"\t",$aa,"\n";
#	print $BB,"\t",$Bb,"\t",$bb,"\n";
#	print  $x2;exit;
	my $line=join("\t", @inf);
	 print C  $line ,"\t",$x2,"\t",$p2,"\t",$Bb,"\t",$Aa,"\n";
##filter
	if ($Isselfedseed==1)
    {
         if($Bb-$Aa>1.001)
         {##80%
        	print E "$inf[0]\t$inf[1]\t",$x2,"\t",$p2,"\t",$Bb,"\t",$Aa,"\n";
         }
	     else
         {
             print D  $line,"\n";
         }
     }
     else
     {
         if($x2>$cuoff)  #....3 99%
         {
              print E "$inf[0]\t$inf[1]\t",$x2,"\t",$p2,"\t",$Bb,"\t",$Aa,"\n";
         }
         else
         {
              print D  $line,"\n";
         }
     }
}
close B;
close C;
close D;
close E;
