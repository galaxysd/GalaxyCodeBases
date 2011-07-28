#!/usr/bin/perl -w

# http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=200987345&clade=mammal&org=Human&db=hg19&hgta_group=allTables&hgta_track=snp132Common&hgta_table=0&hgta_regionType=genome&position=chr21%3A33%2C031%2C597-33%2C041%2C570&hgta_outputType=primaryTable&hgta_outFileName=
# Choose Group=All Tables; table=refGene.
# ’≈—©√∑ 17:42:37
#/ifs1/STT_POP/USER/zhangxm/Project/BrownRat/download/SNP/refgene_to_gff.pl

use strict;
unless(@ARGV)
{
	die "$0 <input file> <output file> <sourceID>\n";
}
open IN,"$ARGV[0]" or die "can't open the input file:$!";
open OUT,">$ARGV[1]" or die "can't write the output:$!";
while(<IN>)
{
	chomp;
	if(/^#/)
	{
		next;
	}
	my @line=split /\t/,$_;
	if($line[6] == $line[7])
	{
		next;
	}
  if($line[-2] ne 'cmpl')
  {
      next;
  }
  if($line[-3] ne 'cmpl')
  {
      next;
  }
  $line[4]++;
  $line[6]++;
	print OUT "$line[2]\t$ARGV[2]\tmRNA\t$line[4]\t$line[5]\t.\t$line[3]\t.\tID=$line[1]; name=$line[12];\n";
	my @exon_start=split /,/,$line[9];
  foreach my $mm(0..$#exon_start)
  {
      $exon_start[$mm]++;
  }
	my @exon_end=split /,/,$line[10];
	my @phase=split /,/,$line[-1];
	foreach my $i(0 .. $#exon_start)
	{
		if($exon_end[$i] < $line[6])
		{
			if($line[3] eq "+")
			{
				print OUT "$line[2]\t$ARGV[2]\t5-UTR\t$exon_start[$i]\t$exon_end[$i]\t.\t$line[3]\t.\tParent=$line[1];\n";
			}
			else
			{
				print OUT "$line[2]\t$ARGV[2]\t3-UTR\t$exon_start[$i]\t$exon_end[$i]\t.\t$line[3]\t.\tParent=$line[1];\n";
			}
			if($exon_start[$i+1])
			{
				my $tem_start=$exon_end[$i]+1;
				my $tem_end=$exon_start[$i+1]-1;
				print OUT "$line[2]\t$ARGV[2]\tintron\t$tem_start\t$tem_end\t.\t$line[3]\t.\tParent=$line[1];\n";
			}
		}
		elsif(($exon_end[$i] >= $line[6])&&($exon_end[$i] <=$line[7]))
		{
			if($exon_start[$i] < $line[6])
			{
				my $tem_start=$exon_start[$i];
				my $tem_end=$line[6]-1;
				if($line[3] eq "+")
				{
					print OUT "$line[2]\t$ARGV[2]\t5-UTR\t$tem_start\t$tem_end\t.\t$line[3]\t.\tParent=$line[1];\n";
				}
				else
				{
					print OUT "$line[2]\t$ARGV[2]\t3-UTR\t$tem_start\t$tem_end\t.\t$line[3]\t.\tParent=$line[1];\n";
				}
				print OUT "$line[2]\t$ARGV[2]\tCDS\t$line[6]\t$exon_end[$i]\t.\t$line[3]\t$phase[$i]\tParent=$line[1];\n";
				if($exon_start[$i+1])
				{
					my $tem_start1=$exon_end[$i]+1;
					my $tem_end1=$exon_start[$i+1]-1;
					print OUT "$line[2]\t$ARGV[2]\tintron\t$tem_start1\t$tem_end1\t.\t$line[3]\t.\tParent=$line[1];\n";
				}
			}
			else
			{
				print OUT "$line[2]\t$ARGV[2]\tCDS\t$exon_start[$i]\t$exon_end[$i]\t.\t$line[3]\t$phase[$i]\tParent=$line[1];\n";
				if($exon_start[$i+1])
				{
					my $tem_start1=$exon_end[$i]+1;
					my $tem_end1=$exon_start[$i+1]-1;
					print OUT "$line[2]\t$ARGV[2]\tintron\t$tem_start1\t$tem_end1\t.\t$line[3]\t.\tParent=$line[1];\n";
				}
			}
		}
		else
		{
			if(($exon_start[$i] <= $line[7])&&($exon_start[$i] >= $line[6]))
			{
				print OUT "$line[2]\t$ARGV[2]\tCDS\t$exon_start[$i]\t$line[7]\t.\t$line[3]\t$phase[$i]\tParent=$line[1];\n";
				my $tem_start=$line[7]+1;
				my $tem_end=$exon_end[$i];
				if($line[3] eq "+")
				{
					print OUT "$line[2]\t$ARGV[2]\t3-UTR\t$tem_start\t$tem_end\t.\t$line[3]\t.\tParent=$line[1];\n";
				}
				else 
				{
					print OUT "$line[2]\t$ARGV[2]\t5-UTR\t$tem_start\t$tem_end\t.\t$line[3]\t.\tParent=$line[1];\n";
				}
				if($exon_start[$i+1])
				{
					my $tem_start=$exon_end[$i]+1;
					my $tem_end=$exon_start[$i+1]-1;
					print OUT "$line[2]\t$ARGV[2]\tintron\t$tem_start\t$tem_end\t.\t$line[3]\t.\tParent=$line[1];\n";
				}
			}
			elsif($exon_start[$i] < $line[6])
			{
				my $tem_start1=$exon_start[$i];
				my $tem_end1=$line[6]-1;
				if($line[3] eq "+")
				{
					print OUT "$line[2]\t$ARGV[2]\t5-UTR\t$tem_start1\t$tem_end1\t.\t$line[3]\t.\tParent=$line[1];\n";
				}
				else
				{
					print OUT "$line[2]\t$ARGV[2]\t3-UTR\t$tem_start1\t$tem_end1\t.\t$line[3]\t.\tParent=$line[1];\n";
				}
				print OUT "$line[2]\t$ARGV[2]\tCDS\t$line[6]\t$line[7]\t.\t$line[3]\t$phase[$i]\tParent=$line[1];\n";
				my $tem_start2=$line[7]+1;
				my $tem_end2=$exon_end[$i];
				if($line[3] eq "+")
				{
					print OUT "$line[2]\t$ARGV[2]\t3-UTR\t$tem_start2\t$tem_end2\t.\t$line[3]\t.\tParent=$line[1];\n";
				}
				else
				{
					print OUT "$line[2]\t$ARGV[2]\t5-UTR\t$tem_start2\t$tem_end2\t.\t$line[3]\t.\tParent=$line[1];\n";
				}
				if($exon_start[$i+1])
				{
					my $tem_start=$exon_end[$i]+1;
					my $tem_end=$exon_start[$i+1]-1;
					print OUT "$line[2]\t$ARGV[2]\tintron\t$tem_start\t$tem_end\t.\t$line[3]\t.\tParent=$line[1];\n";
				}
			}
			else
			{
				if($line[3] eq "+")
				{
					print OUT "$line[2]\t$ARGV[2]\t3-UTR\t$exon_start[$i]\t$exon_end[$i]\t.\t$line[3]\t.\tParent=$line[1];\n";
				}
				else
				{
					print OUT "$line[2]\t$ARGV[2]\t5-UTR\t$exon_start[$i]\t$exon_end[$i]\t.\t$line[3]\t.\tParent=$line[1];\n";
				}
				if($exon_start[$i+1])
				{
					my $tem_start=$exon_end[$i]+1;
					my $tem_end=$exon_start[$i+1]-1;
					print OUT "$line[2]\t$ARGV[2]\tintron\t$tem_start\t$tem_end\t.\t$line[3]\t.\tParent=$line[1];\n";
				}
			}
		}
	}
}
close OUT;
close IN;
