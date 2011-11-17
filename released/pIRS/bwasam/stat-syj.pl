#!/usr/bin/perl -w

# /ifs1/ST_ASMB/USER/shiyj/work/simulate/stat/stat.pl

use strict;
use Getopt::Long;
##get options from command line into variables and set default values
my ($read_len, $in_file, $snp_file, $stat_out, $help);
GetOptions(
	"l:i"=>\$read_len,
	"i:s"=>\$in_file,
	"s:s"=>\$snp_file,
	"o:s"=>\$stat_out,
	"help"=>\$help
);

if($help || !defined $read_len || !defined $in_file)
{
	print "Name

stat.pl  --stat matrix of (ref_base,cycle,read_base,quality) through PE soap file

Usage

 perl stat.pl [option]
 
 	--l\t<num>	read length
 	--i\t<str>	soap.pair
 	--s\t<str>	SNP file [default:none]
 	--o\t<str>	output stat file
 	--h\t\thelp
";
	exit;
}

$stat_out = $in_file . "\.stat" if(!$stat_out);

my $Q_shift=64;

my %SNP;
if($snp_file)
{
	open IN,$ARGV[2] || die "$!";
	while(<IN>)
	{
		chomp;
		my @line=split;
		$SNP{$line[0]}{$line[1]}=1;
	}
	close IN;
}
my $read1_len=$read_len;
my $read2_len=$read1_len;

if($in_file=~/\.gz$/)
{
	open IN,"gzip -dc $in_file|" || die "$!";
}
else
{
	open IN,$in_file || die "$!";
}

my($read1,$read2,@line1,@line2);
#my(%match_stat,%mis_stat);
my(@match,@mis,@quality,$mis_len,$ss,$i,$j);
my %stat;

while($read1=<IN>)
{
	chomp $read1;
	@line1 = split /\s+/,$read1;
	next if($line1[4] ne "a");
	$read2 = <IN>;
	chomp $read2;
	@line2 = split /\s+/,$read2;
	if($line2[4] ne "b")
	{
		print STDERR "The soap result is not pair\n$read1\n$read2\n";
		next;
	} 
	next if($line1[10]=~/S/ || $line2[10]=~/S/);#滤掉部分比对上的
	next if($line1[3]>1 || $line2[3]>1);#滤掉非唯一比对的

	@match = split /[ATCG]+/,$line1[-1];
	@mis = split /\d+/,$line1[-1];
	@quality = split //,$line1[2];
	my $snp;
	if($line1[-1] =~ /^\d/)
	{
		shift @mis;
		if($line1[6] eq "+")
		{
			$ss = 0;
			$snp = 0;
			my @seq = split //,$line1[1];
			for($i=0;$i<@mis;$i++)
			{
				for($j=0;$j<$match[$i];$j++)
				{
					if(!$SNP{$line1[7]}{$line1[8]+$snp})
					{
#						$match_stat{$quality[$ss]}[$ss]++;
						$stat{$seq[$ss]}[$ss]{$seq[$ss]}{$quality[$ss]}++;
					}
					$ss++;
					$snp++;
				}
				$mis_len = length($mis[$i]);
				my @tmp = split //,$mis[$i];
				for($j=0;$j<$mis_len;$j++)
				{
					if(!$SNP{$line1[7]}{$line1[8]+$snp})
					{
#						$mis_stat{$quality[$ss]}[$ss]++;
						$stat{$tmp[$j]}[$ss]{$seq[$ss]}{$quality[$ss]}++;
					}
					$ss++;
					$snp++;
				}
			}
			if($ss<$read1_len)
			{
				for($j=0;$j<$match[$i];$j++)
				{
					if(!$SNP{$line1[7]}{$line1[8]+$snp})
					{
#						$match_stat{$quality[$ss]}[$ss]++;
						$stat{$seq[$ss]}[$ss]{$seq[$ss]}{$quality[$ss]}++;
					}
					$ss++;
					$snp++;
				}
			}
		}
		else
		{
			$ss = $read1_len-1;
			$snp = 0;
			$line1[1]=~tr/ATCG/TAGC/;
			my @seq = split //,$line1[1];
			for($i=0;$i<@mis;$i++)
			{
				for($j=0;$j<$match[$i];$j++)
				{
					if(!$SNP{$line1[7]}{$line1[8]+$snp})
					{
						$stat{$seq[$read1_len-1-$ss]}[$ss]{$seq[$read1_len-1-$ss]}{$quality[$read1_len-1-$ss]}++;
					}
					$ss--;
					$snp++;
				}
				$mis_len = length($mis[$i]);
				$mis[$i] =~ tr/ATCG/TAGC/;
				my @tmp = split //,$mis[$i];
				for($j=0;$j<$mis_len;$j++)
				{
					if(!$SNP{$line1[7]}{$line1[8]+$snp})
					{
						$stat{$tmp[$j]}[$ss]{$seq[$read1_len-1-$ss]}{$quality[$read1_len-1-$ss]}++;
					}
					$ss--;
					$snp++;
				}
			}
			if($ss>=0)
			{
				for($j=0;$j<$match[$i];$j++)
				{
					if(!$SNP{$line1[7]}{$line1[8]+$snp})
					{
						$stat{$seq[$read1_len-1-$ss]}[$ss]{$seq[$read1_len-1-$ss]}{$quality[$read1_len-1-$ss]}++;
					}
					$ss--;
					$snp++;
				}
			}
		}
	}
	else
	{
		die "not all mismatch information start with number";
	}
	@match = split /[ATCG]+/,$line2[-1];
	@mis = split /\d+/,$line2[-1];
	@quality = split //,$line2[2];
	if($line2[-1] =~ /^\d/)
	{
		shift @mis;
		if($line2[6] eq "+")
		{
			$ss = $read1_len;
			$snp = 0;
			my @seq = split//,$line2[1];
			for($i=0;$i<@mis;$i++)
			{
				for($j=0;$j<$match[$i];$j++)
				{
					if(!$SNP{$line2[7]}{$line2[8]+$snp})
					{
						$stat{$seq[$ss-$read1_len]}[$ss]{$seq[$ss-$read1_len]}{$quality[$ss-$read1_len]}++;
					}
					$ss++;
					$snp++;
				}
				$mis_len = length($mis[$i]);
				my @tmp = split //,$mis[$i];
				for($j=0;$j<$mis_len;$j++)
				{
					if(!$SNP{$line2[7]}{$line2[8]+$snp})
					{
						$stat{$tmp[$j]}[$ss]{$seq[$ss-$read1_len]}{$quality[$ss-$read1_len]}++;
					}
					$ss++;
					$snp++;
				}
			}
			if($ss<$read1_len+$read2_len)
			{
				for($j=0;$j<$match[$i];$j++)
				{
					if(!$SNP{$line2[7]}{$line2[8]+$snp})
					{
						$stat{$seq[$ss-$read1_len]}[$ss]{$seq[$ss-$read1_len]}{$quality[$ss-$read1_len]}++;
					}
					$ss++;
					$snp++;
				}
			}
		}
		else
		{
			$ss = $read1_len+$read2_len-1;
			$snp = 0;
			$line2[1] =~ tr/ATCG/TAGC/;
			my @seq = split //,$line2[1];
			for($i=0;$i<@mis;$i++)
			{
				for($j=0;$j<$match[$i];$j++)
				{
					if(!$SNP{$line2[7]}{$line2[8]+$snp})
					{
						$stat{$seq[$read1_len+$read2_len-1-$ss]}[$ss]{$seq[$read1_len+$read2_len-1-$ss]}{$quality[$read1_len+$read2_len-1-$ss]}++;
					}
					$ss--;
					$snp++;
				}
				$mis_len = length($mis[$i]);
				$mis[$i] =~ tr/ATCG/TAGC/;
				my @tmp = split //,$mis[$i];
				for($j=0;$j<$mis_len;$j++)
				{
					if(!$SNP{$line2[7]}{$line2[8]+$snp})
					{
						$stat{$tmp[$j]}[$ss]{$seq[$read1_len+$read2_len-1-$ss]}{$quality[$read1_len+$read2_len-1-$ss]}++;
					}
					$ss--;
					$snp++;
				}
			}
			if($ss>=$read1_len)
			{
				for($j=0;$j<$match[$i];$j++)
				{
					if(!$SNP{$line2[7]}{$line2[8]+$snp})
					{
						$stat{$seq[$read1_len+$read2_len-1-$ss]}[$ss]{$seq[$read1_len+$read2_len-1-$ss]}{$quality[$read1_len+$read2_len-1-$ss]}++;
					}
					$ss--;
					$snp++;
				}
			}
		}
	}
	else
	{
		die "not all mismatch information start with number";
	}
}
close IN;

open OUT,">$stat_out" || die "$!";
my(@qual,@cycle,@ref,@base);
for($i=0;$i<41;$i++)
{
	$qual[$i] = chr($i+$Q_shift);
}
for($i=0;$i<$read1_len+$read2_len;$i++)
{
	$cycle[$i]=$i;
}
$ref[0]="A";
$ref[1]="C";
$ref[2]="G";
$ref[3]="T";
$base[0]="A";
$base[1]="C";
$base[2]="G";
$base[3]="T";

my($r,$c,$b,$q);
print OUT "#ref_base\tcycle\t";
foreach $b (@base)
{
	for($i=0;$i<41;$i++)
	{
		print OUT "$b$i\t";
	}
}
print OUT "\n";
foreach $r (@ref)
{
	foreach $c (@cycle)
	{
		my $c_tmp = $c + 1;
		print OUT "$r\t$c_tmp\t";
		foreach $b (@base)
		{
			foreach $q (@qual)
			{
				$stat{$r}[$c]{$b}{$q}=0 if(!$stat{$r}[$c]{$b}{$q});
				print OUT "$stat{$r}[$c]{$b}{$q}\t";
			}
		}
		print OUT "\n";
	}
}

close OUT;
