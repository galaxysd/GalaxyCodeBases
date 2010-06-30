#!/usr/bin/perl -w
use strict;
unless(@ARGV)
{
	die "<fasta file><depth file><window><output file>\n";
}
open AGCT,"$ARGV[0]"
     or die "can't open the fasta file:$!\n";
open DEPTH,"$ARGV[1]"
     or die "can't open depth file:$!\n";
open OUT,">$ARGV[3]"
     or die "the result can't write:$!\n";

my $win=$ARGV[2];
my $total_num=0;
my $sum=$win;
my ($GC_num,$N_num,$gc_rate,$avg_depth,$seat,$chrom);
my $sum_depth=0;
my (@gcs,@depth,@Ns,@line);
my (%hash_seq,%hash_num,%hash);


while(<AGCT>)
{
	chomp;
	if(/(^>.*$)/&&$total_num<$sum)
	{
		if($total_num==0)
		{
			$chrom=$1;
			$seat=$chrom."+"."$sum";
			next;
		}

		$seat=$chrom."+"."$sum";
		$gc_rate=int(100*$GC_num/($win-$N_num));
	  @{$hash_seq{$seat}}=($N_num,$gc_rate);
		$total_num=0;
		$sum=$win;
		$GC_num=0;
		$N_num=0;
	  $chrom=$1;
		next;
	}
	my $len=length($_);
	$total_num+=$len;
	if($total_num<$sum)
	{
		$GC_num+=$_=~s/(G|C)/$1/gi;
		$N_num+=$_=~s/(N)/$1/gi;
	}
	else
	{
		$GC_num+=substr($_, 0, $sum-$total_num+$len)=~s/(G|C)/$1/gi;
		$N_num+=substr($_, 0, $sum-$total_num+$len)=~s/(N)/$1/gi;
		if($win-$N_num)
		{
			$gc_rate=int(100*$GC_num/($win-$N_num));
		}
		$seat=$chrom."+"."$sum";
		@{$hash_seq{$seat}}=($N_num,$gc_rate);
    $GC_num=substr($_,$sum-$total_num+$len,$total_num-$sum+1)=~s/(G|C)/$1/gi;
    $N_num=substr($_,$sum-$total_num+$len,$total_num-$sum+1)=~s/(N)/$1/gi;
    $sum+=$win;
	}
}
if($total_num<$sum)
{
	$seat=$chrom."+"."$sum";
	$gc_rate=int(100*$GC_num/($win-$N_num));
	@{$hash_seq{$seat}}=($N_num,$gc_rate);
}
#foreach my $key(keys %hash_seq)
#{
#	print "@{$hash_seq{$key}}\n";
#}

while(<DEPTH>)
{
	chomp;
	if(/(^>.*$)/&&$total_num<$sum)
	{
		if($total_num==0)
		{
			$chrom=$1;
			$sum_depth=0;
			$seat=$chrom."+"."$sum";
			next;
		}

		$seat=$chrom."+"."$sum";
    if((exists $hash_seq{$seat})&&($total_num-$sum+$win-${$hash_seq{$seat}}[0])!=0)
    {
    	push @{$hash{${$hash_seq{$seat}}[1]}},$sum_depth/($total_num-$sum+$win-${$hash_seq{$seat}}[0]);
    }
		$total_num=0;
		$sum_depth=0;
		$sum=$win;
	  $chrom=$1;
		next;
	}
	my ($len,$temp_len);      #$len表示$_这一行的数字的总和，$temp_len表示$_这一行数字的个数
	while(/(\d+)/g)
	{
		$len+=$1;
		$temp_len+=1;
	}
	if($total_num<$sum)
	{
		$sum_depth+=$len;
		$total_num+=$temp_len;
	}
	else
	{
		@line=split;
		for(my $i=0;$i<$sum-$total_num+$temp_len;$i++)
		{
	    $sum_depth+=$line[$i];
	  }
		$seat=$chrom."+"."$sum";
    if((exists $hash_seq{$seat})&&($win-${$hash_seq{$seat}}[0])!=0)
    {
    	push @{$hash{${$hash_seq{$seat}}[1]}},$sum_depth/($win-${$hash_seq{$seat}}[0]);
    }
    $sum_depth=0;
    for(my $i=$sum-$total_num+$len;$i<$temp_len;$i++)
    {
    	$sum_depth+=$line[$i];
    }
    $sum+=$win;
	}
}
if($total_num<$sum)
{
	$seat=$chrom."+"."$sum";
  if((exists $hash_seq{$seat})&&($total_num-$sum+$win-${$hash_seq{$seat}}[0]+1)!=0)
  {
  	push @{$hash{${$hash_seq{$seat}}[1]}},$sum_depth/($total_num-$sum+$win-${$hash_seq{$seat}}[0]+1);
  }
}
#foreach my $key(keys %hash)
#{
#	print "@{$hash{$key}}\n";
#}

#输出结果
my @gc_key;           #@gc_key存放gc值
foreach my $key(sort keys %hash)
{
	push @gc_key,$key;
}
print OUT "@gc_key\n";
my $mark=0;
while(1)
{
	for(my $i=0;$i<=$#gc_key;$i++)
	{
		if(@{$hash{$gc_key[$i]}}) #此处改动了原来为$hash{$gc_key[$i]}||$hash{$gc_key[$i]}==0
		{
			my $tem=shift @{$hash{$gc_key[$i]}} ;
		  $mark+=@{$hash{$gc_key[$i]}};
			print OUT $tem,"\t";
		}
		else
		{
			print OUT "NA","\t";
		}
	}
	if($mark>0)
	{
		print OUT "\n";
	  $mark=0;
  }
  else
  {
  	last;
  }
}

close DEPTH;
close OUT;
close AGCT;

__END__
/panfs/GAG/zhangxm/stat/depth.pl
张雪梅

考虑了N区
