#!/usr/bin/perl

=head1 Name



=head1 Description



=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage
  perl draw_sperm.pl  <chr_len_file> <chr_gap_file> <chr_centromere_file> <sperms_halotype_dir>  <sperms_depth_dir>  <sperm_purity_dir>
  --window <int>   window size, i.e number of genomic bases in each window, default=500000
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
use lib "/gpfsdata/Analysis/fanwei/bin/";
use SVG;

##get options from command line into variables and set default values
my ($Window_size, $Verbose,$Help);
GetOptions(
	"window:i"=>\$Window_size,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Window_size ||= 500000;
die `pod2text $0` if (@ARGV == 0 || $Help);

my $chr_len_file = shift;
my $chr_gap_file = shift;
my $chr_centromere_file = shift;
my $phased_snp_dir = shift;
my $sperms_depth_dir = shift;
my $sperms_purity_dir = shift;
my $cross_over_file = shift;
my $sexType_AvgDepth_file = shift;

my $Chrom_num = 24;
my $Sperm_num = 99;

my %ChrLen;
my %ChrGap;  ##存储非gap区，临时改变
my %Centro;
my %SpColor; ##store the sperm haplotype information
my %SpCover; ##store the sperm coverage depth information
my %SpPurity; ##store the sperm purity information
my %CrossOver;
my %CrossOverStatus;

my %StatHaplo;
my %StatSNPdens;
my %StatCovDepth;
my %StatPurity;
my %StatGap;

my %SexType;
my %AvgDepth;

read_chrlen($chr_len_file, \%ChrLen);
read_chrgap($chr_gap_file, \%ChrGap);
read_centromere($chr_centromere_file, \%Centro);
read_sperm_haplo();
read_sperm_depth();
read_sperm_purity();
read_cross_over($cross_over_file);
read_sex_depth($sexType_AvgDepth_file);

my $left_edge_pxi = 30;
my $right_edge_pxi = 30;
my $top_edge_pxi = 30;
my $bottom_edge_pxi = 30;

my $Max_win_num = int($ChrLen{"chr1"} / $Window_size) + 1; ##chr1 is the largest chromosome

my $chr_x_shift = 50;
my $win_y_shift = 1;
my $figure_width = $left_edge_pxi + $right_edge_pxi + $Chrom_num*$chr_x_shift;
my $figure_height = $top_edge_pxi + $bottom_edge_pxi + $win_y_shift*$Max_win_num;
my $chr_width = $chr_x_shift * 0.15;
my $Max_win_snp_num = $Window_size * 0.0005;  ##0.0005 (het-SNP-rate) * 0.8 (call-out rate)
my $Max_win_depth =  2;  ##averge*3; Each sperms sperately.

##draw figures for each sperm
for (my $i=0; $i<$Sperm_num; $i++) {
	my $k = $i + 1;
	my $sperm_id = ($i < 9) ? "S0$k" : "S$k";
	print STDERR "draw sperm ".$sperm_id."\n";
	
	my $sperm_avg_depth = $AvgDepth{$sperm_id};
	my $sperm_SNP_ratio = 0.0;
	my $DiSNP_total = 1390544;
	my $sperm_avg_purity = 0.0;
	my $sperm_avg_purity_divide = 0;

	my $svg = SVG->new('width',$figure_width,'height',$figure_height);
	$svg->rect('x',0, 'y',0,'width',$figure_width,'height',$figure_height,'fill',"white");
		
	for (my $c=1; $c<=$Chrom_num; $c++) {
		my $chr = "chr$c";
		$chr =~ s/23/X/;
		$chr =~ s/24/Y/;

		my $x = $left_edge_pxi + ($c-1)*$chr_x_shift;
		my $y;
		my $text_x = $x;
		my $text_y = $top_edge_pxi-5;
		#my $transform_format = "rotate(-90 $text_x,$text_y)";
		$svg->text('x',$text_x ,'y',$text_y,'fill','black', 'font-family','Arial','font-size',16, '-cdata',$chr);
		my $Window_num = int($ChrLen{$chr} / $Window_size) + 1;
	
		
		##draw the automatic identified crossover sites, as background
		my $cross_p = $CrossOver{$sperm_id}{$chr};
		my $status_p = $CrossOverStatus{$sperm_id}{$chr};
		foreach my $pos (@$cross_p) {
			my $cy = $top_edge_pxi + int($pos/$Window_size)*$win_y_shift;
			my $cx = $x+$chr_width/2;
			my $r = $chr_width*1.2;
			my $status = shift @$status_p;
			if ($status == 1) {
				$svg->line('x1',$cx-$r,'y1',$cy-$r,'x2',$cx+$r,'y2',$cy+$r,'stroke','purple','stroke-width',3);
				$svg->line('x1',$cx-$r,'y1',$cy+$r,'x2',$cx+$r,'y2',$cy-$r,'stroke','purple','stroke-width',3);
			}
			elsif($status == 0) {
				#$svg->line('x1',$cx-$r,'y1',$cy-$r,'x2',$cx+$r,'y2',$cy+$r,'stroke','purple','stroke-width',1);
				#$svg->line('x1',$cx-$r,'y1',$cy+$r,'x2',$cx+$r,'y2',$cy-$r,'stroke','purple','stroke-width',1);
			}
		}

		##draw sperm haplotype
		##use yellow color to represent un-phased region, either caused by low-coverage or conflict links
		my $x = $left_edge_pxi + ($c-1)*$chr_x_shift;
		my $Sp_Chr_p = $SpColor{$chr}[$i];
		my @SNPdensity;
		for (my $win_id = 0; $win_id < $Window_num; $win_id++) {
			my $p = $Sp_Chr_p->[$win_id];
			$p->[0] ||= 0;  ##if it is empty, set default value 0
			$p->[1] ||= 0;  ##if it is empty, set default value 0
			$sperm_SNP_ratio += $p->[0] + $p->[1];
			my $color;
			my $ratio;
			if ($p->[0] >= $p->[1] && $p->[0] >= 2) {
				$color = ($p->[1] / $p->[0] < 0.34) ? "green" : "yellow";
				$ratio = $p->[1] / $p->[0];
			}
			elsif($p->[1] >= $p->[0] && $p->[1] >= 2) {
				$color = ($p->[0] / $p->[1] < 0.34) ? "red" : "yellow";
				$ratio = $p->[0] / $p->[1];
			}
			else{
				$color = "yellow";
				$ratio =  $p->[1] / $p->[0] if($p->[0] > $p->[1]); 
				$ratio =  $p->[0] / $p->[1] if($p->[1] > $p->[0]); 
				$ratio = 1 if($p->[0] == $p->[1]); 
			}
			if ($chr eq "chrX") {
				$color = ($SexType{$sperm_id} eq "X") ? "red" : "white";
			}
			if ($chr eq "chrY") {
				$color = ($SexType{$sperm_id} eq "Y") ? "green" : "white";
			}
			$y = $top_edge_pxi + $win_id*$win_y_shift;
			$svg->rect('x',$x, 'y',$y,'width',$chr_width,'height',$win_y_shift,'stroke',"none", 'fill',$color);
			
			my $ratio_percent = int($ratio * 100);
			$StatHaplo{$sperm_id}{$chr}{$win_id} = "$color($p->[0],$p->[1],$ratio_percent%)";
			
			push (@SNPdensity, $p->[0] + $p->[1]);
		}
		
		##draw gap region in the chromosome
		my $gap_p = $ChrGap{$chr};
		for (my $win_id = 0; $win_id < $Window_num; $win_id++) {
			my $gap_ratio = 1 - $gap_p->[$win_id] / $Window_size;
			my $is_gap = 0;
			if($gap_ratio > 0.9){  ##超过90%为gap的认为整个window为gap,无信号
				$y = $top_edge_pxi + $win_id*$win_y_shift;
				$svg->rect('x',$x, 'y',$y,'width',$chr_width,'height',$win_y_shift,'stroke',"none", "fill", "white");
				$is_gap = 1;
			}
			my $ratio_percent = int($gap_ratio * 100);
			$StatGap{$sperm_id}{$chr}{$win_id} = "$is_gap,$ratio_percent%";
		}

		##draw chromosome frame
		$y = $top_edge_pxi;
		$svg->rect('x',$x, 'y',$y,'width',$chr_width,'height',$Window_num*$win_y_shift,'stroke',"black", "fill", "none");


		##draw centromere (draw on the top of figure)
		my ($centro_start,$centro_end) = ($1,$2) if($Centro{$chr} =~ /(\d+),(\d+)/);
		my $centro_pos = int(($centro_start + $centro_end) / 2);
		$y = $top_edge_pxi + int($centro_pos/$Window_size)*$win_y_shift;
		$svg->circle('cx',$x+$chr_width/2,'cy',$y,'r',$chr_width*0.6,'stroke','black','fill','none','stroke-width',3);

		##draw SNP density curve $Max_win_snp_num
		my $x_start = $left_edge_pxi + ($c-1+0.2)*$chr_x_shift;
		my $x_end = $left_edge_pxi + ($c-1+0.8)*$chr_x_shift;
		my $x_width = $x_end - $x_start;
		$StatSNPdens{$sperm_id}{$chr}{0} = $SNPdensity[0];
		for (my $j=1; $j<@SNPdensity; $j++) {
			my $rate1 = $SNPdensity[$j-1] / $Max_win_snp_num;
			$rate1 = 1 if($rate1 > 1);
			my $x1 = $x_start + $rate1*$x_width;
			my $y1 = $top_edge_pxi + ($j-1)*$win_y_shift;
			my $rate2 = $SNPdensity[$j] / $Max_win_snp_num;
			$rate2 = 1 if($rate2 > 1);
			my $x2 = $x_start + $rate2*$x_width;
			my $y2 = $top_edge_pxi + $j*$win_y_shift;
			$svg->line('x1',$x1,'y1',$y1,'x2',$x2,'y2',$y2,'stroke','orange','stroke-width',1);
			
			$StatSNPdens{$sperm_id}{$chr}{$j} = $SNPdensity[$j];
		}
		
		##draw coverage depth curve
		$Max_win_depth = int($AvgDepth{$sperm_id}*3*10)/10;
		my $depth_p = $SpCover{$chr}[$i];
		$StatCovDepth{$sperm_id}{$chr}{0} = ($gap_p->[0] == 0) ? 0 : int($depth_p->[0] / $gap_p->[0] * 100) / 100;
		for (my $j=1; $j<$Window_num; $j++) {
			my $rate1 = ($gap_p->[$j-1] == 0) ? 0 : ($depth_p->[$j-1] / $gap_p->[$j-1] / $Max_win_depth);
			$rate1 = 1 if($rate1 > 1);
			my $x1 = $x_start + $rate1*$x_width;
			my $y1 = $top_edge_pxi + ($j-1)*$win_y_shift;
			my $rate2 = ($gap_p->[$j] == 0) ? 0 : ($depth_p->[$j] / $gap_p->[$j] / $Max_win_depth);
			$rate2 = 1 if($rate2 > 1);
			my $x2 = $x_start + $rate2*$x_width;
			my $y2 = $top_edge_pxi + $j*$win_y_shift;
			$svg->line('x1',$x1,'y1',$y1,'x2',$x2,'y2',$y2,'stroke','blue','stroke-width',1);
			
			$StatCovDepth{$sperm_id}{$chr}{$j} = ($gap_p->[$j] == 0) ? 0: int($depth_p->[$j] / $gap_p->[$j]*100)/100;
		}
		
		##draw purity curve
		my $purity_p = $SpPurity{$chr}[$i];
		for (my $j=1; $j<$Window_num; $j++) {
			my $pre_p = $purity_p->[$j-1];
			my $cur_p = $purity_p->[$j];
			$pre_p->[0] ||= 0;
			$pre_p->[1] ||= 0;
			$cur_p->[0] ||= 0;
			$cur_p->[1] ||= 0;
			my ($rate1,$rate2);
			my $pre_total = $pre_p->[0] + $pre_p->[1];
			my $cur_total = $cur_p->[0] + $cur_p->[1];
			if ($pre_total <= 0) {
				$rate1 = 0;
			}else{
				$rate1 = ($pre_p->[0] > $pre_p->[1]) ? ($pre_p->[0] / $pre_total) : ($pre_p->[1] / $pre_total);
			}
			if ($cur_total <= 0) {
				$rate2 = 0;
			}else{
				$rate2 = ($cur_p->[0] > $cur_p->[1]) ? ($cur_p->[0] / $cur_total) : ($cur_p->[1] / $cur_total);
				$sperm_avg_purity += $rate2;
				$sperm_avg_purity_divide ++;
			}
			my $x1 = $x_start + $rate1*$x_width;
			my $y1 = $top_edge_pxi + ($j-1)*$win_y_shift;
			my $x2 = $x_start + $rate2*$x_width;
			my $y2 = $top_edge_pxi + $j*$win_y_shift;
			
			$svg->line('x1',$x1,'y1',$y1,'x2',$x2,'y2',$y2,'stroke','cyan','stroke-width',1);
			
			if ($j == 1) {
				my $rate_percent = int($rate1*100);
				my $num1 = int($pre_p->[0] * 100)/100;
				my $num2 = int($pre_p->[1] * 100)/100;
				$StatPurity{$sperm_id}{$chr}{0} = "$num1,$num2,$rate_percent%";
			}
			my $rate_percent = int($rate2*100);
			my $num1 = int($cur_p->[0] * 100)/100;
			my $num2 = int($cur_p->[1] * 100)/100;
			$StatPurity{$sperm_id}{$chr}{$j} = "$num1,$num2,$rate_percent%";
		
		}
		
		##draw small ruler scale for each chromosome properties
		my $x1 = $x_start + 0*$x_width;
		my $y1 = $top_edge_pxi + 0*$win_y_shift;
		my $x2 = $x_start + 1*$x_width;
		my $y2 = $top_edge_pxi + $Window_num*$win_y_shift;
		#$svg->line('x1',$x1,'y1',$y1,'x2',$x2,'y2',$y1,'stroke','black','stroke-width',1);
		#$svg->line('x1',$x1,'y1',$y2,'x2',$x2,'y2',$y2,'stroke','black','stroke-width',1);
		##画个小尺子，分成10个格即可，并非真正表示坐标
		#plot_ruler("svg",$svg,"Y",$y1, "X_start",$x1,"X_end",$x2,"bp_start",1,"bp_end",1000000,"scaletype","Mb","scaletypepos","left","scalestart","force","rulerstyle",1,"bigscalesize",1000000);
		plot_ruler("svg",$svg,"Y",$y2, "X_start",$x1,"X_end",$x2,"bp_start",1,"bp_end",1000000,"scaletype","Mb","scaletypepos","left","scalestart","force","rulerstyle",2,"bigscalesize",1000000);

		$svg->text('x',$figure_width/4,'y',$figure_height-$bottom_edge_pxi-10,'fill','black', 'font-family','Arial','font-size',24, '-cdata', "Sperm ID: $sperm_id");
	}
	
	##draw the legends
	my $x = $figure_width* 0.7;
	my $y = $figure_height / 2;
	my $skip_y = $win_y_shift*$Max_win_num / 2 / 9;
	my $x_skip = $chr_x_shift;
	my $x_text = $x + $x_skip + 10;
	my $legend_font_size = 16;
	$svg->rect('x',$x, 'y',$y,'width',$x_skip,'height',$chr_width,'stroke',"black", "fill", "green");
	$svg->text('x',$x_text,'y',$y+$chr_width,'fill','black', 'font-family','Arial','font-size',$legend_font_size, '-cdata', "Paternal haplotype");
	$y += $skip_y;
	$svg->rect('x',$x, 'y',$y,'width',$x_skip,'height',$chr_width,'stroke',"black", "fill", "red");
	$svg->text('x',$x_text,'y',$y+$chr_width,'fill','black', 'font-family','Arial','font-size',$legend_font_size, '-cdata', "Maternal haplotype");
	$y += $skip_y;
	$svg->rect('x',$x, 'y',$y,'width',$x_skip,'height',$chr_width,'stroke',"black", "fill", "yellow");
	$svg->text('x',$x_text,'y',$y+$chr_width,'fill','black', 'font-family','Arial','font-size',$legend_font_size, '-cdata', "Unresolved haplotype");
	$y += $skip_y;
	$svg->rect('x',$x, 'y',$y,'width',$x_skip,'height',$chr_width,'stroke',"black", "fill", "white");
	$svg->text('x',$x_text,'y',$y+$chr_width,'fill','black', 'font-family','Arial','font-size',$legend_font_size, '-cdata', "Unassembled gaps in chromosome");
	$y += $skip_y;
	$svg->circle('cx',$x+$x_skip/2,'cy',$y,'r',$chr_width*0.6,'stroke','black','fill','none','stroke-width',3);
	$svg->text('x',$x_text,'y',$y+$chr_width,'fill','black', 'font-family','Arial','font-size',$legend_font_size, '-cdata', "Centromere region");
	$y += $skip_y;
	
	my $cx = $x+$x_skip/2;
	my $cy = $y;
	my $r = $chr_width*1.2;
	$svg->line('x1',$cx-$r,'y1',$cy-$r,'x2',$cx+$r,'y2',$cy+$r,'stroke','purple','stroke-width',3);
	$svg->line('x1',$cx-$r,'y1',$cy+$r,'x2',$cx+$r,'y2',$cy-$r,'stroke','purple','stroke-width',3);
	$svg->text('x',$x_text,'y',$y+$chr_width,'fill','black', 'font-family','Arial','font-size',$legend_font_size, '-cdata', "Manual curated crossover");
	$y += $skip_y;

	$svg->line('x1',$x,'y1',$y,'x2',$x+$x_skip,'y2',$y,'stroke','blue','stroke-width',1);
	$svg->text('x',$x_text,'y',$y+$chr_width,'fill','black', 'font-family','Arial','font-size',$legend_font_size, '-cdata', "Coverage depth (Max: $Max_win_depth)");
	$y += $skip_y;
	$svg->line('x1',$x,'y1',$y,'x2',$x+$x_skip,'y2',$y,'stroke','orange','stroke-width',1);
	$svg->text('x',$x_text,'y',$y+$chr_width,'fill','black', 'font-family','Arial','font-size',$legend_font_size, '-cdata', "Density of called SNP (Max: 5E-4)");
	$y += $skip_y;
	$svg->line('x1',$x,'y1',$y,'x2',$x+$x_skip,'y2',$y,'stroke','cyan','stroke-width',1);
	$svg->text('x',$x_text,'y',$y+$chr_width,'fill','black', 'font-family','Arial','font-size',$legend_font_size, '-cdata', "Purity of covered reads (Max: 100%)");
	
	##draw overall stat. numbers
	$sperm_avg_depth = int($sperm_avg_depth*100)/100;
	$sperm_SNP_ratio =  int($sperm_SNP_ratio / $DiSNP_total * 100)/100;
	$sperm_avg_purity = int($sperm_avg_purity / $sperm_avg_purity_divide *100)/100;
	my $win_size = ($Window_size / 1000)." kb";
	$x = $figure_width* 0.45;
	$y = $figure_height * 0.7;
	$skip_y = $win_y_shift*$Max_win_num * 0.3 / 4;
	$legend_font_size = 16;
	
	$svg->text('x',$x,'y',$y,'fill','black', 'font-family','Arial','font-size',$legend_font_size, '-cdata', "Avg. coverage depth: $sperm_avg_depth");
	$y += $skip_y;
	$svg->text('x',$x,'y',$y,'fill','black', 'font-family','Arial','font-size',$legend_font_size, '-cdata', "Ratio of called SNP: $sperm_SNP_ratio");
	$y += $skip_y;
	$svg->text('x',$x,'y',$y,'fill','black', 'font-family','Arial','font-size',$legend_font_size, '-cdata', "Overall purity: $sperm_avg_purity");
	$y += $skip_y;
	$svg->text('x',$x,'y',$y,'fill','black', 'font-family','Arial','font-size',$legend_font_size, '-cdata', "Stat. window size: $win_size");
	$y += $skip_y;

	my $svg_file = "$sperm_id.allchrs.haplotype.snp.depth.purity.svg";
	open OUT, ">$svg_file" || die "fail create $svg_file\n";
	print OUT $svg->xmlify();
	close OUT;
	`svg2xxx.pl  -memory 20G  $svg_file`;
}


##output the result numbers in the svg figure　
##make statistics for the whole genome and chromosomes
for (my $i=0; $i<$Sperm_num; $i++) {
	my $k = $i + 1;
	my $sperm_id = ($i < 9) ? "S0$k" : "S$k";
	
	open OUT, ">$sperm_id.allchrs.haplotype.snp.depth.purity.res";
	print OUT "#SpID\tChr\tWinStart\tWinSize\tIsGap\tHaplo\tSNPnum\tCovDep\tPurity\n";
	for (my $c=1; $c<=$Chrom_num; $c++) {
		my $chr = "chr$c";
		$chr =~ s/23/X/;
		$chr =~ s/24/Y/;
		my $Window_num = int($ChrLen{$chr} / $Window_size) + 1;
		for (my $j=0; $j<$Window_num; $j++) {
			my $window_start = $j*$Window_size + 1;
			print OUT "$sperm_id\t$chr\t$window_start\t$Window_size\t".$StatGap{$sperm_id}{$chr}{$j}."\t".$StatHaplo{$sperm_id}{$chr}{$j}."\t".$StatSNPdens{$sperm_id}{$chr}{$j}."\t".$StatCovDepth{$sperm_id}{$chr}{$j}."\t".$StatPurity{$sperm_id}{$chr}{$j}."\n";
		}
	}
	close OUT;
}


##organize the so-many generated result files
mkdir("sperms-haplotype-figures-png-win$Window_size");
mkdir("sperms-haplotype-figures-svg-win$Window_size");
mkdir("sperms-haplotype-figures-res-win$Window_size");
`mv *.allchrs.haplotype.snp.depth.purity.png sperms-haplotype-figures-png-win$Window_size`;
`mv *.allchrs.haplotype.snp.depth.purity.svg sperms-haplotype-figures-svg-win$Window_size`;
`mv *.allchrs.haplotype.snp.depth.purity.res sperms-haplotype-figures-res-win$Window_size`;

##make final summary of the coverage depth and purity for chromosomes and genome
my %SumCovDepth;
my %SumPurity;
my %WinNum;
my $Total_win_num = 0;

for (my $c=1; $c<=$Chrom_num; $c++) {
	my $chr = "chr$c";
	$chr =~ s/23/X/;
	$chr =~ s/24/Y/;
	my $Window_num = int($ChrLen{$chr} / $Window_size) + 1;
	for (my $j=0; $j<$Window_num; $j++) {
		if ($StatGap{"S01"}{$chr}{$j} =~ /^0,/) {  ##只统计非gap区的
			$WinNum{$chr} ++;
			$WinNum{"Total"} ++;
		}
	}
}

for (my $i=0; $i<$Sperm_num; $i++) {
	my $k = $i + 1;
	my $sperm_id = ($i < 9) ? "S0$k" : "S$k";
	
	for (my $c=1; $c<=$Chrom_num; $c++) {
		my $chr = "chr$c";
		$chr =~ s/23/X/;
		$chr =~ s/24/Y/;
		my $Window_num = int($ChrLen{$chr} / $Window_size) + 1;
		for (my $j=0; $j<$Window_num; $j++) {
			if ($StatGap{$sperm_id}{$chr}{$j} =~ /^0,/) {  ##只统计非gap区的
				my $cov_depth = $StatCovDepth{$sperm_id}{$chr}{$j};
				my $purity_val = $1/100 if($StatPurity{$sperm_id}{$chr}{$j} =~ /(\d+).$/);
				$SumCovDepth{$sperm_id}{$chr} += $cov_depth;
				$SumCovDepth{$sperm_id}{"Total"} += $cov_depth;
				$SumPurity{$sperm_id}{$chr} += $purity_val;
				$SumPurity{$sperm_id}{"Total"} += $purity_val;
			}
		}
	}
	
}


my $head = "#SpID";
for (my $c=1; $c<=$Chrom_num; $c++) {
	my $chr = "chr$c";
	$chr =~ s/23/X/;
	$chr =~ s/24/Y/;
	$head .= "\t$chr";
}
$head .= "\tTotal";
my $out_depth = $head;
my $out_purity = $head;

for (my $i=0; $i<$Sperm_num; $i++) {
	my $k = $i + 1;
	my $sperm_id = ($i < 9) ? "S0$k" : "S$k";
	$out_depth .= "\n$sperm_id";
	$out_purity .= "\n$sperm_id";
	for (my $c=1; $c<=$Chrom_num; $c++) {
		my $chr = "chr$c";
		$chr =~ s/23/X/;
		$chr =~ s/24/Y/;
		$out_depth .= ($WinNum{$chr} == 0) ? "\t0" : "\t".($SumCovDepth{$sperm_id}{$chr} / $WinNum{$chr});
		$out_purity .= ($WinNum{$chr} == 0) ? "\t0" : "\t".($SumPurity{$sperm_id}{$chr} / $WinNum{$chr});
	}
	$out_depth .= ($WinNum{"Total"} == 0) ? "\t0" : "\t".($SumCovDepth{$sperm_id}{"Total"} / $WinNum{"Total"});
	$out_purity .= ($WinNum{"Total"} == 0) ? "\t0" : "\t".($SumPurity{$sperm_id}{"Total"} / $WinNum{"Total"});
}

open OUT, ">all.sperms.cover.depth.stat" || die "fail";
print OUT $out_depth;
close OUT;
open OUT, ">all.sperms.cover.purity.stat" || die "fail";
print OUT $out_purity;
close OUT;

####################################################
################### Sub Routines ###################
####################################################

sub read_chrlen{
	my $file = shift;
	my $ChrLen_p = shift;

	open IN, $file || die "fail $file";
	while (<IN>) {
		next if(/^\#/);
		chomp;
		my ($chr,$len) = ($1,$2) if(/^(\S+)\s+(\d+)/);
		$ChrLen_p->{$chr} = $len;
	}
	close IN;
}

sub read_chrgap{
	my $file = shift;
	my $ChrGap_p = shift;

	open IN, $file || die "fail $file";
	while (<IN>) {
		next if(/^\#/);
		chomp;
		my @t = split /\s+/;
		my $chr = $t[0];
		my $pos = $t[1];
		my $bases_num = $t[4];
		my $win_id = int($pos / $Window_size);

		$ChrGap_p->{$chr}[$win_id] += $bases_num;
	}
	close IN;
}


##23      chr1    121535434       124535434       1270    N       3000000 centromere      no
sub read_centromere{
	my $file = shift;
	my $Centro_p = shift;

	open IN, $file || die "fail $file";
	while (<IN>) {
		next if(/^\#/);
		chomp;
		my @t = split /\s+/;
		my $chr = $t[1];
		my $start = $t[2];
		my $end = $t[3];
		$Centro_p->{$chr} = "$start,$end";
	}
	close IN;
}

##reading data into memory %SpColor, read the phased *.filter.haplo.more.FaMo (~80% of all het SNP)
sub read_sperm_haplo{
	
	my @HapSpermFiles = glob("$phased_snp_dir/chr*.diploid.sperms.snp.strict.filter.haplo.more.FaMo");

	foreach my $haplotype_sperm_file (@HapSpermFiles) {
		
		next if($haplotype_sperm_file =~ /chrX/ || $haplotype_sperm_file =~ /chrY/);
		
		##reading SNP data for each sperms
		print STDERR "reading file ".$haplotype_sperm_file."\n";
		
		##calculate the window data
		##四种颜色：红和绿为亲本之一，黄色为数据足够但混乱区；白色为缺乏数据区。
		my $Chr_name;

		open IN, $haplotype_sperm_file || die "fail $haplotype_sperm_file";
		while (<IN>) {
			next if(/^\#/);
			chomp;
			my @t = split /\t/;
			$Chr_name = $t[0];
			my $pos = $t[1];
			my $hap1 = $t[2];
			my $hap2 = $t[3];
			my $SpNum = $t[4];
			my $win_id = int($pos/$Window_size);
			for (my $i=0; $i<$Sperm_num; $i++) {
				my $sperm_base = $t[$i+5];
				if ($sperm_base eq $hap1) {
					$SpColor{$Chr_name}[$i][$win_id][0] ++;  ##0 represent green
				}
				elsif($sperm_base eq $hap2) {
					$SpColor{$Chr_name}[$i][$win_id][1] ++;  ##1 represent red
				}
			}
		}
		close IN;
		
	}

}

##reading the covered number of reads for the phased SNP, to check the purity of sperm
sub read_sperm_purity{
	
	my @Files = glob("$sperms_purity_dir/chr*.diploid.sperms.phased_snp.cover.depth");

	foreach my $sperm_purity_file (@Files) {
		
		next if($sperm_purity_file =~ /chrX/ || $sperm_purity_file =~ /chrY/);
		
		##reading SNP data for each sperms
		print STDERR "reading file ".$sperm_purity_file."\n";
		
		##calculate the window data
		##四种颜色：红和绿为亲本之一，黄色为数据足够但混乱区；白色为缺乏数据区。
		my $Chr_name;

		open IN, $sperm_purity_file || die "fail $sperm_purity_file";
		while (<IN>) {
			next if(/^\#/);
			chomp;
			my @t = split /\t/;
			$Chr_name = $t[0];
			my $pos = $t[1];
			my $hap1 = $t[2];
			my $hap2 = $t[3];
			my $SpNum = $t[4];
			my $win_id = int($pos/$Window_size);
			for (my $i=0; $i<$Sperm_num; $i++) {
				my $sperm_cover = $t[$i+5];
				my ($base1,$depth1,$base2,$depth2) = ($1,$2,$3,$4) if($sperm_cover =~ /(\w),(\d+);(\w),(\d+)/); 
				my $depth_total = $depth1 + $depth2;
				next if($base1 eq "N" || $base2 eq "N" || $depth_total <= 0);
				my $rate1 = $depth1 / $depth_total;
				my $rate2 = $depth2 / $depth_total;
				if ($base1 eq $hap1) {
					$SpPurity{$Chr_name}[$i][$win_id][0] += $rate1;
					$SpPurity{$Chr_name}[$i][$win_id][1] += $rate2;
				}
				if ($base2 eq $hap1) {
					$SpPurity{$Chr_name}[$i][$win_id][0] += $rate2;
					$SpPurity{$Chr_name}[$i][$win_id][1] += $rate1;
				}
			}
		}
		close IN;
		
	}

}


##reading data into memory %SpCover
sub read_sperm_depth{
	
	my @DepSpermFiles = glob("$sperms_depth_dir/chr*.allsperms.depth.stat");
	
	foreach my $sperm_depth_file (@DepSpermFiles) {
		
		next if($sperm_depth_file =~ /chrM/);
		
		##reading SNP data for each sperms
		print STDERR "reading file ".$sperm_depth_file."\n";
		
		open IN, $sperm_depth_file || die "fail $sperm_depth_file";
		while (<IN>) {
			next if(/^\#/);
			chomp;
			my @t = split /\t/;
			my $Chr_name = shift @t;
			my $pos = shift @t;
			my $win_id = int($pos / $Window_size);
			shift @t;
			for (my $i=0; $i<@t; $i++) {
				$SpCover{$Chr_name}[$i][$win_id] += $t[$i];  
			}

		}
		close IN;
	
	}

}
#ruler stype 3, 刻度朝上
##get configurations from a hash
##the following fields must be specified: $rulcfg{svg},$rulcfg{Y},$rulcfg{X_start},$rulcfg{X_end},$rulcfg{bp_start},$rulcfg{bp_end}
##the following fields can be auto-assigned: $rulcfg{scaletype}, $rulcfg{scaletypepos}, $rulcfg{scalestart}, $rulcfg{rulerstyle},$rulcfg{font_family}, $rulcfg{font_size}
##usage example: plot_ruler("svg",$svg,"Y",$backbone1_height - $y_margin/3, "X_start",$synteny_Xstart,"X_end",$synteny_Xstart+$Seq1_len/$synteny_length*$synteny_width,"bp_start",$Seq1_start,"bp_end",$Seq1_end,"scaletype","bp","scaletypepos","right","scalestart","force","rulerstyle",3,"bigscalesize",1000000);
sub plot_ruler{
	my %rulcfg = @_;
	$rulcfg{scaletype} ||= "bp";
	$rulcfg{scaletypepos} ||= "left";
	$rulcfg{scalestart} ||= "auto";
	$rulcfg{rulerstyle} ||= "1";
	$rulcfg{font_family} ||= "ArialNarrow";
	$rulcfg{font_size} ||= 32;
	
	
	my $scale_size = 6;
	my $divid = 50;
	my $unit_start;

	my $bp_len = $rulcfg{bp_end} - $rulcfg{bp_start};
	my $X_len = $rulcfg{X_end} - $rulcfg{X_start};

	##caculate the length of smallest unit 
	my ($str,$str1,$str2,$unit);
	$str = $bp_len / $divid;
	$str = sprintf("%e",$str);
	if ($str=~/([\d\.]+)e([+-\d]+)/) {
		$str1 = $1;
		$str2 = $2;
	}
	$str1 = int ( $str1 + 0.5 );
	$unit = $str1 * 10 ** $str2;
	$unit = $rulcfg{bigscalesize}/10 if(defined $rulcfg{bigscalesize});

	my $g = $rulcfg{svg}->group('id'=>times().rand());
	
	## draw the main axis
	$g->line('x1',$rulcfg{X_start},'y1',$rulcfg{Y},'x2',$rulcfg{X_end},'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1);		
	return if($rulcfg{bp_end} - $rulcfg{bp_start}  == 0);
	
#	##draw ruler mark text at the specified postion(left or right of the ruler)
#	if ($rulcfg{scaletypepos} eq "left") {
#		$g->text('x',$rulcfg{X_start}-textWidth($rulcfg{font_family},$rulcfg{font_size},$rulcfg{scaletype})-6,'y',$rulcfg{Y},'-cdata',$rulcfg{scaletype},"font-family",$rulcfg{font_family},"font-size",$rulcfg{font_size},"fill",'#000000');
#	}
#	if ($rulcfg{scaletypepos} eq "right") {
#		$g->text('x',$rulcfg{X_end} + 6,'y',$rulcfg{Y},'-cdata',$rulcfg{scaletype},"font-family",$rulcfg{font_family},"font-size",$rulcfg{font_size},"fill",'#000000');
#	}

	##decide unit start
	if ($rulcfg{scalestart} eq "auto") {
		$unit_start = $rulcfg{bp_start} + ($unit - $rulcfg{bp_start} % $unit);
	}
	if ($rulcfg{scalestart} eq "force") {
		$unit_start = int($rulcfg{bp_start} / 10 + 0.5) * 10; ##四舍五入，零碎取整
	}
	

	## draw small scale lines
	for (my $i=$unit_start; $i<=$rulcfg{bp_end}; $i+=$unit) {
		my $X = $rulcfg{X_start} + ($i - $rulcfg{bp_start}) / $bp_len * $X_len;
		$g->line('x1',$X,'y1',$rulcfg{Y} - $scale_size/2,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '1');
		$g->line('x1',$X,'y1',$rulcfg{Y} + $scale_size/2,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '2');
		$g->line('x1',$X,'y1',$rulcfg{Y} - $scale_size/2,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '3');
	}
	## draw big scale lines and text scales
	for (my $i=$unit_start; $i<=$rulcfg{bp_end}; $i+=$unit*10) {
		my $X = $rulcfg{X_start} + ($i - $rulcfg{bp_start}) / $bp_len * $X_len;
		$g->line('x1',$X,'y1',$rulcfg{Y} - $scale_size,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '1');
		$g->line('x1',$X,'y1',$rulcfg{Y} + $scale_size,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1)  if ($rulcfg{rulerstyle} eq '2');
		$g->line('x1',$X,'y1',$rulcfg{Y} - $scale_size,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '3');
		my $disp_scale_text = $i / $rulcfg{bigscalesize} if(defined $rulcfg{bigscalesize});
		#$rulcfg{svg}->text('x',$X - textWidth($rulcfg{font_family},$rulcfg{font_size},$disp_scale_text) / 2,'y',$rulcfg{Y}+textHeight($rulcfg{font_size}),'fill','#000000','-cdata',$disp_scale_text,'font-size',$rulcfg{font_size}, 'font-family',$rulcfg{font_family}) if ($rulcfg{rulerstyle} eq '1');
		#$rulcfg{svg}->text('x',$X - textWidth($rulcfg{font_family},$rulcfg{font_size},$disp_scale_text) / 2,'y',$rulcfg{Y}+$scale_size+textHeight($rulcfg{font_size})+2,'fill','#000000','-cdata',$disp_scale_text,'font-size',$rulcfg{font_size}, 'font-family',$rulcfg{font_family})  if ($rulcfg{rulerstyle} eq '2');
		#$rulcfg{svg}->text('x',$X - textWidth($rulcfg{font_family},$rulcfg{font_size},$disp_scale_text) / 2,'y',$rulcfg{Y}-$scale_size-textHeight($rulcfg{font_size})+6,'fill','#000000','-cdata',$disp_scale_text,'font-size',$rulcfg{font_size}, 'font-family',$rulcfg{font_family}) if ($rulcfg{rulerstyle} eq '3');
	}

}


sub read_cross_over(){
	my $file = shift;

	open IN, $file || die "fail ";
	while (<IN>) {
		next if(/^\#/);
		chomp;
		my @t = split /\s+/;
		my $sperm_id = $t[0];
		my $chr = $t[1];

		next if($chr eq "chrX");
		my $pos = $t[2];
		my $curate_status = $t[6];
		$curate_status =~ s/\s//g;
		if ($curate_status != 1 && $curate_status != 0) {
			print STDERR "Error in the crossover file, $curate_status\n";
			exit();
		}
		if ($curate_status !~ /1/ && $curate_status !~ /0/) {
			print STDERR "Error in the crossover file, $curate_status\n";
			exit();
		}

		push @{$CrossOver{$sperm_id}{$chr}}, $pos;
		push @{$CrossOverStatus{$sperm_id}{$chr}}, $curate_status;
	}
	close IN;
}


sub read_sex_depth{
	my $file = shift;
	open IN, $file || die "fail ";
	while (<IN>) {
		next if(/^\#/ || /^SpId/);
		chomp;
		my @t = split /\s+/;
		my $sperm_id = $t[0];
		my $sextype = $t[-1];
		my $avgdepth = $t[-2];
		$SexType{$sperm_id} = $sextype;
		$AvgDepth{$sperm_id} = $avgdepth;
	}
	close IN;
}
