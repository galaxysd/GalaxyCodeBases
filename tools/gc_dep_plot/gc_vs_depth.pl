#!/usr/bin/perl -w
#use strict;
use Time::HiRes qw ( gettimeofday tv_interval );
use Getopt::Long;
my ($depth_file,$fa_file,$out_file,$win_size,$help);
GetOptions(
	"depth:s"=>\$depth_file,
	"fa:s"=>\$fa_file,
	"output:s"=>\$out_file,
	"win:i"=>\$win_size,
	"help"=>\$help
);
if(!$depth_file || !$fa_file || !$out_file || !$win_size || $help)
{
	print "\t -depth: the depth file\n";
	print "\t -fa: the .fa file\n";
	print "\t -output: the output file\n";
	print "\t -win: the window size\n";
	print "\t -help: if you need the help";
	die;
}
#set time
my $str_time=[gettimeofday];
#getDepth
print "Begin to count average of depth\n";
my %depth=();
my $cnt=0;
my $sum=0;
my $chr;
open(IN,"<$depth_file")or die "can not open $depth_file :$!\n";
while(<IN>)
{
	chomp;
	if($_!~/^>/)
	{
		while($_=~/(\d+)/g)
		{
			my $num=$1;
			if($num==65535){$num=0;}
			$sum+=$num;
			$cnt++;
			if($cnt==$win_size)
			{
				my $aver=$sum/$cnt;
				push @{$chr},$aver;
				$sum=0;
				$cnt=0;
			}
		}
	}
	else
	{
		$chr=$_;
		$depth{$chr}=$chr;
		print "$chr\n";
	}
}
close IN;
close OUT;
print "Finish the depth counting\n";

#countGC
$chr=();
$cnt=0;
$sum=0;
print "Begin to count the gc\n";
open(IN,"<$fa_file")or die "can not open $fa_file:$!";
while(<IN>)
{
	chomp;
	if($_!~/^>/)
	{
		my $num=s/[gc]/Y/ig;
		$sum+=$num;
		$cnt+=length($_);
		if($cnt==$win_size)
		{
			my $aver=int($sum/$cnt*100);
			push @{$chr},$aver;
			$sum=0;
			$cnt=0;
		}
		elsif($cnt>$win_size)
		{
			my $dis=$cnt-$win_size;
			my @chr=split(//,$_);
			my $len=@chr;
			my $count=0;
			for(my $i=$len-$dis;$i<$len;$i++)
			{
				if($chr[$i]=~/[gc]/i)
				{
					$count++;
				}
			}
			my $aver=int(($sum-$count)/$win_size*100);
                	push @{$chr},$aver;
                	$sum=$count;
                	$cnt=$dis;
		}
		else {next;}
	}
	else
	{
		$chr="gc$_";
		print "$_\n";
	}
}
close IN;

foreach(keys %depth)
{
	my $chrname=$_;
	foreach(@{"gc$chrname"})
	{
		my $dep=shift @{$chrname};
		push @{$_},$dep;
	}
}

my $str=0;
my $end=80;
my $na=0;
my $line=$end-$str;
my $max_depth=0;
my $cnt_depth=0;
open(OUT,">$out_file")or die "$!\n";
for(my $i=$str;$i<=$end;$i++)
{
	print OUT "$i\t";
}
print OUT "\n";
while($na<=$line)
{
	$na=0;
	for(my $j=$str;$j<=$end;$j++)
	{
		my $dep=shift @{$j};
		if(!defined $dep){$dep="NA";$na++;}
		else{$max_depth+=$dep;$cnt_depth++;}
		print OUT "$dep\t";
	}
	print OUT "\n";
}
close OUT;
#R plot
$max_depth/=$cnt_depth/2;
$max_depth=(int($max_depth/50)+1)*50;
my $area=0;
foreach(1..$max_depth/50)
{
	my $ev_depth=$_*50;
	$area.=",$ev_depth";
}
open(OUT,">$out_file.R")or die;
print OUT 'a=read.delim("',"$out_file",'");',"\n";
print OUT 'pdf("',"$out_file.pdf",'",width=8,height=6,paper="special");',"\n";
print OUT "boxplot(a,ylim=c(0,$max_depth),xlim=c(0,90),axes=FALSE,outline=FALSE,xlab='GC%',ylab='Depth');\n";
print OUT "axis(1,at=c(0,10,20,30,40,50,60,70,80),labels=c(0,10,20,30,40,50,60,70,80),las=0);\n";
print OUT "axis(2,at=c($area),labels=c($area),las=2);\nq();\nn\n";
close OUT;
`Rscript $out_file.R`;
my $end_time=[gettimeofday];
my $time=tv_interval($str_time,$end_time);
