#!/usr/bin/perl

=head1 Name



=head1 Description



=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage
  
  perl infer_haplotype.pl  <*.sperms.snp.file>
  --snpnum <int>  set the minimum number of identified sperm SNPs, default=20
  --snpscore <int> set the minimum score for each sperm SNP, default=10
  --minor <float>  the number of each identified sperm should be >= $num_spsnp/$Minor_ratio_times, default=3
  --verbose   output running progress information to screen  
  --help      output help information to screen  

=head1 Exmple



=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;
use File::Path;  ## function " mkpath" and "rmtree" deal with directory

##get options from command line into variables and set default values
my ($Minor_ratio_times,$Snpnum_Cutoff, $Snpscore_cutoff, $Link_ratio, $Verbose,$Help);
GetOptions(
	"snpnum:i"=>\$Snpnum_Cutoff,
	"snpscore:i"=>\$Snpscore_cutoff,
	"minor:f"=>\$Minor_ratio_times,
	"link:f"=>\$Link_ratio,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Snpnum_Cutoff ||= 20;
$Snpscore_cutoff ||= 10;
$Minor_ratio_times ||= 3;
$Link_ratio ||= 0.2;
die `pod2text $0` if (@ARGV == 0 || $Help);

my $sperms_snp_file = shift;

my @DiSNP;
my @SpermSNP;
my $CHR;
my @Pos;
my $Header;

open IN, $sperms_snp_file || die "fail $sperms_snp_file";
while (<IN>) {
	if(/^\#/){
		chomp;
		$Header = $_;
		next;
	}
	my @t = split /\t/;
	my $chr = shift @t;
	$CHR = $chr;
	my $pos = shift @t;
	my $disnp = shift @t;
	my ($base1, $base2, $disnp_depth, $disnp_score) = ($1,$2,$3,$4) if($disnp =~ /(\w)(\w)\),(\d+),(\d+)$/);
	my $num_spsnp = 0;
	my @genotypes;
	foreach  (@t) {
		my ($base, $score) = ($1,$2) if(/(\w+),\d+,(\d+)$/);
		if ($score >= $Snpscore_cutoff) {
			$num_spsnp ++;
		}else{
			$base = "N";
		}
		push @genotypes, $base;
	}
	next if($num_spsnp < $Snpnum_Cutoff);
	
	my %hash;
	$hash{$base1} = 0;
	$hash{$base2} = 0;
	foreach my $ba (@genotypes) {
		$hash{$ba} ++;
	}
	if ( $hash{$base1} >= $num_spsnp/$Minor_ratio_times && $hash{$base2} >= $num_spsnp/$Minor_ratio_times ) {
		push @Pos, $pos;
		push @DiSNP, [$base1, $base2];
		push @SpermSNP, \@genotypes;
	}

}
close IN;


#get a good start point fisrt before iteration, 找相邻的两个
my $good_start = 0;
for (my $i=1; $i<@DiSNP; $i++) {
	my $pre_dibase1 = $DiSNP[$i - 1][0];
	my $pre_dibase2 = $DiSNP[$i - 1][1];
	my $cur_dibase1 = $DiSNP[$i][0];
	my $cur_dibase2 = $DiSNP[$i][1];
	
	my $pre_spermP = $SpermSNP[$i - 1];
	my $cur_spermP = $SpermSNP[$i];
	my %links;
	for (my $j=0; $j<@$cur_spermP; $j++) {
		$links{$pre_spermP->[$j]}{$cur_spermP->[$j]} ++;
	}
	my $link1 = $links{$pre_dibase1}{$cur_dibase1} || 0;
	my $link2 = $links{$pre_dibase2}{$cur_dibase2} || 0;
	my $link3 = $links{$pre_dibase1}{$cur_dibase2} || 0;
	my $link4 = $links{$pre_dibase2}{$cur_dibase1} || 0;
	my $ratio;
	my $status;
	
	if ($link3 + $link4 >= $link1 + $link2) {
		 $ratio = ($link3 + $link4 > 0) ? (($link1 + $link2) / ($link3 + $link4)) : 1;
		$status = 1;
	}else{
		$ratio = ($link3 + $link4) / ($link1 + $link2);
		$status = -1;
	}
	
	if ($ratio <= $Link_ratio/2) {
		$good_start = $i;
		last;
	}
}


print STDERR "good_start: $good_start in $sperms_snp_file\n";


#print $Header,"link1_num\tlink2_num\tlink3_num\tlink4_num\tlink_ratio\n";
my $last_i = $good_start - 1;

my (@HapPos, @Hap1, @Hap2);
push @HapPos, $Pos[$last_i];
push @Hap1, $DiSNP[$last_i][0];
push @Hap2, $DiSNP[$last_i][1];
my $pre_status = 1;

for (my $i=$good_start; $i<@DiSNP; $i++) {
	my $pre_dibase1 = $DiSNP[$last_i][0];
	my $pre_dibase2 = $DiSNP[$last_i][1];
	my $cur_dibase1 = $DiSNP[$i][0];
	my $cur_dibase2 = $DiSNP[$i][1];
	
	my $pre_spermP = $SpermSNP[$last_i];
	my $cur_spermP = $SpermSNP[$i];
	my %links;
	for (my $j=0; $j<@$cur_spermP; $j++) {
		$links{$pre_spermP->[$j]}{$cur_spermP->[$j]} ++;
	}
	my $link1 = $links{$pre_dibase1}{$cur_dibase1} || 0;
	my $link2 = $links{$pre_dibase2}{$cur_dibase2} || 0;
	my $link3 = $links{$pre_dibase1}{$cur_dibase2} || 0;
	my $link4 = $links{$pre_dibase2}{$cur_dibase1} || 0;
	
	my $ratio;
	my $status;
	
	if ($link3 + $link4 >= $link1 + $link2) {
		 $ratio = ($link3 + $link4 > 0) ? (($link1 + $link2) / ($link3 + $link4)) : 1;
		$status = -1;
	}else{
		$ratio = ($link3 + $link4) / ($link1 + $link2);
		$status = 1;
	}

	

	if ($ratio <= $Link_ratio) {
		push @HapPos, $Pos[$i];
		if ($status == $pre_status) {
			push @Hap1, $DiSNP[$i][0];
			push @Hap2, $DiSNP[$i][1];
			$pre_status = 1;
		}else{
			push @Hap1, $DiSNP[$i][1];
			push @Hap2, $DiSNP[$i][0];
			$pre_status = -1;
		}

		my $line = $CHR."\t".$Pos[$i]."\t".$Hap1[-1]."\t".$Hap2[-1];
		for (my $j=0; $j<@$cur_spermP; $j++) {			
			$line .= "\t".$cur_spermP->[$j];
		}
		my $snp_dist = $i - $last_i;
		$line .= "\t$link1\t$link2\t$link3\t$link4\t$ratio\t$snp_dist\n";
		
		print $line;
		$last_i = $i;
	}
	
}

