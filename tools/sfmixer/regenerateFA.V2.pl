#! /usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
#use File::Basename;

my %opts = ();

GetOptions(\%opts, "faFile:s", "len:i", "insert:i", "outDir:s", "oriChrName:s", "chrname:s", "project:s");

unless (defined $opts{faFile} && $opts{outDir} && $opts{project}) {
	print "\n\tthis program will generate a new fa file based on the original fa file. The scaffolds will be merged into several chromosomes using \"n\" to sperate them\n\tthe output will be the new fa file and a list with which the scaffold will be traced back\n";
	print "\n\t-faFile\t\tthe original fa file\n";
	print "\t-outDir\t\tthe output directory\n";
	print "\t-project\tproject name\n";
	print "\t-len\t\tthe len of each new chromosome after merge [10000000]\n";
	print "\t-insert\t\tthe number of \"n\" insert into the new chromosome between two scaffolds [75]\n";
	print "\t-oriChrName\tthe original name pattern (case insensitive) [Chr]\n";
	print "\t-chrname\tthe new chromosome name pattern [ChrNew]\n";
	print "Example:\n\tperl $0 -faFile /share/tmp/pub/Genome/cucumberBGI/Cucumber.original.fa -outDir /share/raid11/zhanghao/software/preProcess/ -project Cucumber -chrname LG_M\n";
	print "\toutput: /share/raid11/zhanghao/software/preProcess/Cucumber.merge.fa\n\t\t/share/raid11/zhanghao/software/preProcess/Cucumber.merge.list\n";
	print "\tperl $0 -faFile /share/tmp/pub/Genome/SorghumBicolor/sbi1.fasta -outDir /share/raid11/zhanghao/software/preProcess/ -project SorghumBicolor -oriChrName chromosome -chrname chromosome_M\n";
	print "\tAuther: Hao Zhang\tTime: 22:32 21/04/2009\n";
	exit 0;
}

my $fa_file = $opts{faFile};
my $chr_len = $opts{len};
my $insert_size = $opts{insert};
my $out_dir = $opts{outDir};
my $scaffold_tem = $opts{oriChrName};
my $chr_tem = $opts{chrname};
my $project_name = $opts{project};

$chr_len ||= 10000000;
my $offset = $chr_len/1000;
$insert_size ||= 75;
$scaffold_tem ||= "Chr";
$chr_tem ||= "ChrNew";
my $base_per_line = 50;
#print"$scaffold_tem\t$chr_tem\n";

my $chr_num = 1;
my $this_chr_index = 1;
my $merging = 0;
my $current_chr = "";
my $current_scaf = "";
my %hChrInfo = ();
my %hScaName = ();	#{chr} => (sca1, sca2, sca3...)

my $starttime = time();

open FA, $fa_file or die "$!";
#my $file_basename = basename $fa_file;
#$file_basename =~ s/\.fa//;
open MERGEFA, ">$out_dir/$project_name.merge.fa" or die "$!";
my $base_per_line_mark = 1;
while (my $line = <FA>) {
	chomp $line;
	if ($line =~ /^>/) {
		if ($line !~ /$scaffold_tem/i) {
			my $scaffold_name = (split /\s+/, $line)[0];
			$scaffold_name =~ s/^>//;
			if ($this_chr_index > ($chr_len - $offset)) { #finish merging one chromosome, and will start merge a new one.
				$hChrInfo{$chr_num}{$current_scaf} .= "\t$this_chr_index";

				++$chr_num;
				$current_chr = "$chr_tem$chr_num";
				$current_scaf = $scaffold_name;
				$this_chr_index = 1;
				$hChrInfo{$chr_num}{$current_scaf} = $this_chr_index;
				push @{$hScaName{$chr_num}}, $current_scaf;
				print MERGEFA "\n\>$current_chr\n";
			}
			else { #still the previous one
				if ((1 == $this_chr_index) and (1 == $chr_num)) { #the first scaffold
					$current_chr = "$chr_tem$chr_num";
					$current_scaf = $scaffold_name;
					$hChrInfo{$chr_num}{$current_scaf} = $this_chr_index;
					push @{$hScaName{$chr_num}}, $current_scaf;
					print MERGEFA ">$current_chr\n";
				}
				else {
					my $end_index = $this_chr_index - 1;
					$hChrInfo{$chr_num}{$current_scaf} .= "\t$end_index";
					for(my $count = $this_chr_index; ($this_chr_index - $count) < $insert_size; ++$this_chr_index) {
						print MERGEFA "n";
						if (0 == ($this_chr_index % $base_per_line)) {
							print MERGEFA "\n";
						}
					}
					$current_scaf = $scaffold_name;
					$hChrInfo{$chr_num}{$current_scaf} = $this_chr_index;
					push @{$hScaName{$chr_num}}, $current_scaf;
				}
			}
			$merging = 1;
		} #end if
		else {
			$merging = 0;
			print MERGEFA "$line\n";
		}
	}
	else {
		if (0 == $merging) {
			print MERGEFA "$line\n";
			if ($base_per_line_mark) {
				my @count_base_per_line = unpack "(a)*", $line;
				$base_per_line = @count_base_per_line;
				$base_per_line_mark = 0;
			}
		}
		else {
			my @bases = unpack "(a)*", $line;
			foreach my $base (@bases) {
				print MERGEFA "$base";
				if(0 == ($this_chr_index % $base_per_line)) {
					print MERGEFA "\n";
				}
				++$this_chr_index;
			}
		}
	}
} #end while
my $end_index = $this_chr_index - 1;
$hChrInfo{$chr_num}{$current_scaf} .= "\t$end_index";
close MERGEFA;
close FA;

open LIST, ">$out_dir/$project_name.merge.list" or die "$!";
open CHRORDER, ">$out_dir/$project_name.merge.chrorder" or die "$!";
foreach my $chr (sort {$a<=>$b} keys %hChrInfo) {
	print CHRORDER "$chr_tem$chr\n";
	foreach my $scaf (@{$hScaName{$chr}}) {
		print LIST "$chr_tem$chr\t$scaf\t$hChrInfo{$chr}{$scaf}\n";
	}
}
close CHRORDER;
close LIST;

my $endtime = time();
my $processtime = $endtime - $starttime;
print "take $processtime seconds to process $fa_file\n";
