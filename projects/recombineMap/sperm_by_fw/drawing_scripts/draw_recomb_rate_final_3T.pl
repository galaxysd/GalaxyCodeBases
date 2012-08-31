#!/usr/bin/perl

=head1 Name



=head1 Description



=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage
  perl draw_recomb_rate.pl  <chr_len_file> <chr_gap_file> <chr_centromere_file> <sperms_recomb_file>  <HapMap_recomb_file>  <deCODE_recomb_file>
  --window <int>   window size, i.e number of genomic bases in each window, default=500000
  --yshift <int>   the pixels for one window size, default=1
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
my ($Window_size, $win_y_shift, $Verbose,$Help);
GetOptions(
	"window:i"=>\$Window_size,
	"yshift:i"=>\$win_y_shift,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$Window_size ||= 500000;
$win_y_shift ||= 1;
my $Window_size_PAR = $Window_size / 20;
die `pod2text $0` if (@ARGV == 0 || $Help);

print STDERR "Window_size: $Window_size\n";
print STDERR "win_y_shift: $win_y_shift\n";

my $chr_len_file = shift;
my $chr_gap_file = shift;
my $chr_centromere_file = shift;
my $sperms_notused_file = shift;
my $sperms_recomb_file = shift;
my $HapMap_recomb_file = shift;
my $deCODE_recomb_file = shift;

my $Chrom_num = 23;
my $Sperm_num = 0;
my $Max_cM_vs_Mb_ratio = 5;
my %ChrLen;
my %ChrGap;  ##存储非gap区，临时改变
my %Centro;
my %CrossOver;
my %NotUsedSperms;

my %StatRecomb;
my %SpRecombRate;
my %HapMapRate;
my %StatHapMap;
my %deCodeRate;
my %StatDeCode;

my %StatAccum;
my %StatCorre;

read_chrlen($chr_len_file, \%ChrLen);
read_chrgap($chr_gap_file, \%ChrGap);
read_centromere($chr_centromere_file, \%Centro);
read_sperm_list($sperms_notused_file);
read_cross_over($sperms_recomb_file);
read_HapMap_file($HapMap_recomb_file);
read_deCODE_file($deCODE_recomb_file);

$Sperm_num = keys %StatRecomb;

print STDERR "The used sperm number is: $Sperm_num\n";

my $left_edge_pxi = 30;
my $right_edge_pxi = 30;
my $top_edge_pxi = 30;
my $bottom_edge_pxi = 30;

my $Max_win_num = int($ChrLen{"chr1"} / $Window_size) + 1; ##chr1 is the largest chromosome

my $chr_x_shift = 50;

my $figure_width = $left_edge_pxi + $right_edge_pxi + $Chrom_num*$chr_x_shift;
my $figure_height = $top_edge_pxi + $bottom_edge_pxi + $win_y_shift*$Max_win_num;
my $chr_width = $chr_x_shift * 0.15;

##draw figures for recombination rate

my $svg = SVG->new('width',$figure_width,'height',$figure_height);
$svg->rect('x',0, 'y',0,'width',$figure_width,'height',$figure_height,'fill',"white");
	
for (my $c=1; $c<=22; $c++) {
	my $chr = "chr$c";
	my $x = $left_edge_pxi + ($c-1)*$chr_x_shift;
	my $y;

	$y = $top_edge_pxi-8;
	$svg->text('x',$x-5 ,'y',$y,'fill','black', 'font-family','Arial','font-size',16, '-cdata',$chr);
	
	my $Window_num = int($ChrLen{$chr} / $Window_size) + 1;
	
	my $x = $left_edge_pxi + ($c-1)*$chr_x_shift;
	
	##draw chromosome frame
	$y = $top_edge_pxi;
	$svg->rect('x',$x, 'y',$y,'width',$chr_width,'height',$Window_num*$win_y_shift,'stroke',"black", "fill", "none");

	##draw gap region in the chromosome
	my $gap_p = $ChrGap{$chr};
	for (my $win_id = 0; $win_id < $Window_num; $win_id++) {
		my $gap_ratio = 1 - $gap_p->[$win_id] / $Window_size;
		my $is_gap = 0;
		if($gap_ratio > 0.9){  ##超过90%为gap的认为整个window为gap,无信号
			$y = $top_edge_pxi + $win_id*$win_y_shift;
			$svg->rect('x',$x, 'y',$y,'width',$chr_width,'height',$win_y_shift,'stroke',"none", "fill", "black");
			$is_gap = 1;
		}
		my $ratio_percent = int($gap_ratio * 100);
	}

	##draw chromosome frame
	$y = $top_edge_pxi;
	$svg->rect('x',$x, 'y',$y,'width',$chr_width,'height',$Window_num*$win_y_shift,'stroke',"black", "fill", "none");


	##draw centromere (draw on the top of figure)
	my ($centro_start,$centro_end) = ($1,$2) if($Centro{$chr} =~ /(\d+),(\d+)/);
	my $centro_pos = int(($centro_start + $centro_end) / 2);
	$y = $top_edge_pxi + int($centro_pos/$Window_size)*$win_y_shift;
	$svg->circle('cx',$x+$chr_width/2,'cy',$y,'r',$chr_width*0.6,'stroke','red','fill','red','stroke-width',1);
	
	##draw the automatic identified crossover sites
	my $cross_p = $CrossOver{$chr};
	foreach my $pos (@$cross_p) {
		my $cy = $top_edge_pxi + int($pos/$Window_size)*$win_y_shift;
		my $cx = $x+$chr_width/2;
		my $r = $chr_width*0.8;
		$svg->line('x1',$cx-$r,'y1',$cy-$r,'x2',$cx+$r,'y2',$cy+$r,'stroke','purple','stroke-width',1);
		$svg->line('x1',$cx-$r,'y1',$cy+$r,'x2',$cx+$r,'y2',$cy-$r,'stroke','purple','stroke-width',1);
	}

	my $x_start = $left_edge_pxi + ($c-1+0.2)*$chr_x_shift;
	my $x_end = $left_edge_pxi + ($c-1+0.8)*$chr_x_shift;
	my $x_width = $x_end - $x_start;
	
	##draw sperm recombination rate curve, the value is cM/Mb
	my $accum_value = 0.0;
	my $Sp_recomb_p = $SpRecombRate{$chr};
	for (my $j=1; $j<$Window_num; $j++) {
		$Sp_recomb_p->[$j-1] ||= 0;
		$Sp_recomb_p->[$j] ||= 0;
		my $rate1 = ( ($Sp_recomb_p->[($j-2>=0) ? ($j-2) : 0] + $Sp_recomb_p->[$j-1] + $Sp_recomb_p->[$j]) / $Sperm_num * 100) / (3*$Window_size / 1000000) / $Max_cM_vs_Mb_ratio;
		$rate1 = 1 if($rate1 > 1);
		my $x1 = $x_start + $rate1*$x_width;
		my $y1 = $top_edge_pxi + ($j-1)*$win_y_shift;
		my $rate2 = ( ($Sp_recomb_p->[$j-1] + $Sp_recomb_p->[$j] + $Sp_recomb_p->[($j+1<$Window_num) ? ($j+1) : $j]) / $Sperm_num * 100) / (3*$Window_size / 1000000) / $Max_cM_vs_Mb_ratio;
		$rate2 = 1 if($rate2 > 1);
		my $x2 = $x_start + $rate2*$x_width;
		my $y2 = $top_edge_pxi + $j*$win_y_shift;
		$svg->line('x1',$x1,'y1',$y1,'x2',$x2,'y2',$y2,'stroke','blue','stroke-width',1);
		
		$accum_value += $Sp_recomb_p->[$j-1] / $Sperm_num * 100;
		$StatAccum{$chr}{$j*$Window_size}{"sperm"} = $accum_value;

		my $gap_ratio = 1 - $gap_p->[$j] / $Window_size;
		if ($gap_ratio > 0 && ($j*$Window_size > 5000000) && ($j*$Window_size < $ChrLen{$chr} - 5000000) ) {
			push @{$StatCorre{$chr}{"sperm"}} , $Sp_recomb_p->[$j-1] / $Sperm_num * 100;
		}
		
	}

	##draw HapMapII recombination rate curve, the value is cM/Mb
	$accum_value = 0.0;
	my $Hapmap_p = $HapMapRate{$chr};
	for (my $j=1; $j<$Window_num; $j++) {
		$Hapmap_p->[$j-1] ||= 0;
		$Hapmap_p->[$j] ||= 0;
		my $rate1 = ( $Hapmap_p->[($j-2>=0) ? ($j-2) : 0] + $Hapmap_p->[$j-1] + $Hapmap_p->[$j] ) / (3*$Window_size / 1000000) / $Max_cM_vs_Mb_ratio;
		$rate1 = 1 if($rate1 > 1);
		my $x1 = $x_start + $rate1*$x_width;
		my $y1 = $top_edge_pxi + ($j-1)*$win_y_shift;
		my $rate2 = ($Hapmap_p->[$j-1] + $Hapmap_p->[$j] +$Hapmap_p->[($j+1<$Window_num) ? ($j+1) : $j]) / (3*$Window_size / 1000000) / $Max_cM_vs_Mb_ratio;
		$rate2 = 1 if($rate2 > 1);
		my $x2 = $x_start + $rate2*$x_width;
		my $y2 = $top_edge_pxi + $j*$win_y_shift;
		$svg->line('x1',$x1,'y1',$y1,'x2',$x2,'y2',$y2,'stroke','orange','stroke-width',1);
		
		$accum_value += $Hapmap_p->[$j-1];
		$StatAccum{$chr}{$j*$Window_size}{"hapmap"} = $accum_value;
		
		my $gap_ratio = 1 - $gap_p->[$j] / $Window_size;
		if ($gap_ratio > 0 && ($j*$Window_size > 5000000) && ($j*$Window_size < $ChrLen{$chr} - 5000000) ) {
			push @{$StatCorre{$chr}{"hapmap"}} , $Hapmap_p->[$j-1];
		}
		
	}

	##draw deCode male recombination rate curve, the value is cM/Mb
	$accum_value = 0.0;
	my $deCode_p = $deCodeRate{$chr};
	for (my $j=1; $j<$Window_num; $j++) {
		$deCode_p->[$j-1] ||= 0;
		$deCode_p->[$j] ||= 0;
		my $rate1 = ($deCode_p->[($j-2>=0) ? ($j-2) : 0] + $deCode_p->[$j-1] + $deCode_p->[$j] ) / (3*$Window_size / 1000000) / $Max_cM_vs_Mb_ratio;
		$rate1 = 1 if($rate1 > 1);
		my $x1 = $x_start + $rate1*$x_width;
		my $y1 = $top_edge_pxi + ($j-1)*$win_y_shift;
		my $rate2 = ($deCode_p->[$j-1] + $deCode_p->[$j] + $deCode_p->[($j+1<$Window_num) ? ($j+1) : $j]) / (3*$Window_size / 1000000) / $Max_cM_vs_Mb_ratio;
		$rate2 = 1 if($rate2 > 1);
		my $x2 = $x_start + $rate2*$x_width;
		my $y2 = $top_edge_pxi + $j*$win_y_shift;
		$svg->line('x1',$x1,'y1',$y1,'x2',$x2,'y2',$y2,'stroke','cyan','stroke-width',1);
		
		$accum_value += $deCode_p->[$j-1];
		$StatAccum{$chr}{$j*$Window_size}{"decode"} = $accum_value;
		
		my $gap_ratio = 1 - $gap_p->[$j] / $Window_size;
		if ($gap_ratio > 0 && ($j*$Window_size > 5000000) && ($j*$Window_size < $ChrLen{$chr} - 5000000) ) {
			push @{$StatCorre{$chr}{"decode"}}, $deCode_p->[$j-1];
		}
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
	$svg->text('x',$x1-4,'y',$y2+15,'fill','black', 'font-family','Arial','font-size',12, '-cdata', "0");
	$svg->text('x',$x2-4,'y',$y2+15,'fill','black', 'font-family','Arial','font-size',12, '-cdata', "$Max_cM_vs_Mb_ratio");
}

##draw the legends
my $x = $figure_width* 0.15;
my $y = $figure_height * 0.85;
my $skip_y = 30;
my $x_skip = $chr_x_shift;
my $x_text = $x + $x_skip + 10;
my $legend_font_size = 16;

$svg->rect('x',$x, 'y',$y,'width',$x_skip,'height',$chr_width,'stroke',"black", "fill", "black");
$svg->text('x',$x_text,'y',$y+$chr_width,'fill','black', 'font-family','Arial','font-size',$legend_font_size, '-cdata', "Unassembled gaps");
$y += $skip_y;
$svg->rect('x',$x, 'y',$y,'width',$x_skip,'height',$chr_width,'stroke',"black", "fill", "white");
$svg->text('x',$x_text,'y',$y+$chr_width,'fill','black', 'font-family','Arial','font-size',$legend_font_size, '-cdata', "Continuous sequences");

$x = $figure_width* 0.35;
$y = $figure_height * 0.7;
$skip_y = 30;
$x_skip = $chr_x_shift;
$x_text = $x + $x_skip + 10;
$legend_font_size = 16;

$svg->circle('cx',$x+$x_skip/2,'cy',$y,'r',$chr_width*0.6,'stroke','red','fill','red','stroke-width',1);
$svg->text('x',$x_text,'y',$y+$chr_width,'fill','black', 'font-family','Arial','font-size',$legend_font_size, '-cdata', "Centromere");
$y += $skip_y;

my $cx = $x+$x_skip/2;
my $cy = $y;
my $r = $chr_width;
$svg->line('x1',$cx-$r,'y1',$cy-$r,'x2',$cx+$r,'y2',$cy+$r,'stroke','purple','stroke-width',1);
$svg->line('x1',$cx-$r,'y1',$cy+$r,'x2',$cx+$r,'y2',$cy-$r,'stroke','purple','stroke-width',1);
$svg->text('x',$x_text,'y',$y+$chr_width,'fill','black', 'font-family','Arial','font-size',$legend_font_size, '-cdata', "Curated crossover");
$y += $skip_y;

$svg->line('x1',$x,'y1',$y,'x2',$x+$x_skip,'y2',$y,'stroke','blue','stroke-width',1);
$svg->text('x',$x_text,'y',$y+$chr_width,'fill','black', 'font-family','Arial','font-size',$legend_font_size, '-cdata', "Sperm cM/Mb");
$y += $skip_y;
$svg->line('x1',$x,'y1',$y,'x2',$x+$x_skip,'y2',$y,'stroke','orange','stroke-width',1);
$svg->text('x',$x_text,'y',$y+$chr_width,'fill','black', 'font-family','Arial','font-size',$legend_font_size, '-cdata', "HapMap cM/Mb");
$y += $skip_y;
$svg->line('x1',$x,'y1',$y,'x2',$x+$x_skip,'y2',$y,'stroke','cyan','stroke-width',1);
$svg->text('x',$x_text,'y',$y+$chr_width,'fill','black', 'font-family','Arial','font-size',$legend_font_size, '-cdata', "deCODE cM/Mb");

#$svg->text('x',$figure_width/4,'y',$figure_height-$bottom_edge_pxi,'fill','black', 'font-family','Arial','font-size',24, '-cdata', "Recombination events and  hotspots");



##output recombination statistics
my $sperm_stat_file = "recombination_by_sperms.tab";
my $chr_stat_file = "recombination_by_chroms.tab";
my $Recomb_total = 0;
my %StatChr;
open OUT, ">$sperm_stat_file" || die "fail";
my $head = "#SpId";
for (my $i=1; $i<=23; $i++) {
	my $chr = "chr$i";
	$chr =~ s/23/X/g;
	$head .= "\t$chr";
}
$head .= "\tTotal";
print OUT $head."\n";
foreach my $spermId (sort keys %StatRecomb) {
	my $sperm_p = $StatRecomb{$spermId};
	my $line = $spermId;
	my $sperm_total = 0;
	for (my $i=1; $i<=23; $i++) {
		my $chr = "chr$i";
		$chr =~ s/23/X/g;
		$sperm_p->{$chr} ||= 0;
		$StatChr{$chr} += $sperm_p->{$chr};
		$Recomb_total += $sperm_p->{$chr};
		$line .= "\t".$sperm_p->{$chr};
		$sperm_total += $sperm_p->{$chr};
	}
	$line .= "\t".$sperm_total;
	print OUT $line."\n";

}
close OUT;


open OUT, ">$chr_stat_file" || die "fail";
print OUT "#Chr\tLen(Mb)\tAvg_recomb(#)\tMap(cM)\tRate(cM/Mb)\tHapMap(cM)\tRate(cM/Mb)\tdeCODEmale(cM)\tRate(cM/Mb)\n";
my $genome_len = 0;
for (my $i=1; $i<=23; $i++) {
	my $chr = "chr$i";
	$chr =~ s/23/X/g;
	my $chr_len = $ChrLen{$chr};
	
	$genome_len += $chr_len;
	my $chr_len_int = int($chr_len / 1000000); ##Mb
	my $avg_recomb = int($StatChr{$chr} / $Sperm_num * 100) / 100 ;
	my $cM_len = $avg_recomb * 100;  ##Mb
	my $cM_phy_ratio = ($StatChr{$chr} / $Sperm_num * 100) / ($chr_len / 1000000);
	my $HapMap_cM = int($StatHapMap{$chr});
	my $HapMap_cMvsMb_rate =  $HapMap_cM / ($chr_len / 1000000);
	my $deCode_cM = int($StatDeCode{$chr});
	my $deCode_cMvsMb_rate = $deCode_cM / ($chr_len / 1000000);
	my $line = "$chr\t$chr_len_int\t$avg_recomb\t$cM_len\t$cM_phy_ratio\t$HapMap_cM\t$HapMap_cMvsMb_rate\t$deCode_cM\t$deCode_cMvsMb_rate\n";
	print OUT $line;
}
my $chr = "Genome";
my $chr_len = int($genome_len / 1000000); ##Mb
my $avg_recomb = int($Recomb_total / $Sperm_num * 100) / 100 ;
my $cM_len = $avg_recomb * 100;  ##Mb
my $cM_phy_ratio = ($Recomb_total / $Sperm_num * 100) / ($genome_len / 1000000);
my $HapMap_total = int($StatHapMap{"Total"});
my $HapMap_total_rate = $StatHapMap{"Total"} / ($genome_len / 1000000);
my $DeCode_total = int($StatDeCode{"Total"});
my $DeCode_total_rate = $StatDeCode{"Total"} / ($genome_len / 1000000);
my $line = "$chr\t$chr_len\t$avg_recomb\t$cM_len\t$cM_phy_ratio\t$HapMap_total\t$HapMap_total_rate\t$DeCode_total\t$DeCode_total_rate\n";
print OUT $line;
close OUT;




##draw chromosome X PAR region, -150M if >= 3M, the same way as sperm chrX PAR.
my $chr = "chrX";
my $x = $left_edge_pxi + 22*$chr_x_shift;
my $y;

$y = $top_edge_pxi-8;
$svg->text('x',$x-5 ,'y',$y,'fill','black', 'font-family','Arial','font-size',14, '-cdata',"PAR1");

$Max_cM_vs_Mb_ratio *= 20;
my $Window_num = int($ChrLen{$chr} / $Window_size_PAR) + 1;
	
##draw chromosome frame
$y = $top_edge_pxi;
$svg->rect('x',$x, 'y',$y,'width',$chr_width,'height',$Window_num*$win_y_shift,'stroke',"black", "fill", "none");

##draw the middle gap between par1 and par2 (left)


##draw chromosome frame
$y = $top_edge_pxi;
$svg->rect('x',$x, 'y',$y,'width',$chr_width,'height',$Window_num*$win_y_shift,'stroke',"black", "fill", "none");

##draw the automatic identified crossover sites
my $cross_p = $CrossOver{$chr};
foreach my $pos (@$cross_p) {
	my $cy = $top_edge_pxi + int($pos/$Window_size_PAR)*$win_y_shift;
	my $cx = $x+$chr_width/2;
	my $r = $chr_width*0.8;
	$svg->line('x1',$cx-$r,'y1',$cy-$r,'x2',$cx+$r,'y2',$cy+$r,'stroke','purple','stroke-width',1);
	$svg->line('x1',$cx-$r,'y1',$cy+$r,'x2',$cx+$r,'y2',$cy-$r,'stroke','purple','stroke-width',1);
}

my $x_start = $left_edge_pxi + (22+0.2)*$chr_x_shift;
my $x_end = $left_edge_pxi + (22+0.8)*$chr_x_shift;
my $x_width = $x_end - $x_start;

##draw sperm recombination rate curve, the value is cM/Mb
my $accum_value = 0.0;
my $Sp_recomb_p = $SpRecombRate{$chr};
for (my $j=1; $j<$Window_num; $j++) {
	$Sp_recomb_p->[$j-1] ||= 0;
	$Sp_recomb_p->[$j] ||= 0;
	my $rate1 = ( ($Sp_recomb_p->[($j-2>=0) ? ($j-2) : 0] + $Sp_recomb_p->[$j-1] + $Sp_recomb_p->[$j]) / $Sperm_num * 100) / (3*$Window_size_PAR / 1000000) / $Max_cM_vs_Mb_ratio;
	$rate1 = 1 if($rate1 > 1);
	my $x1 = $x_start + $rate1*$x_width;
	my $y1 = $top_edge_pxi + ($j-1)*$win_y_shift;
	my $rate2 = (($Sp_recomb_p->[$j-1] + $Sp_recomb_p->[$j] + $Sp_recomb_p->[($j+1<$Window_num) ? ($j+1) : $j]) / $Sperm_num * 100) / (3*$Window_size_PAR / 1000000) / $Max_cM_vs_Mb_ratio;
	$rate2 = 1 if($rate2 > 1);
	my $x2 = $x_start + $rate2*$x_width;
	my $y2 = $top_edge_pxi + $j*$win_y_shift;
	$svg->line('x1',$x1,'y1',$y1,'x2',$x2,'y2',$y2,'stroke','blue','stroke-width',1);
	
	$accum_value += $Sp_recomb_p->[$j-1] / $Sperm_num * 100;
	$StatAccum{$chr}{$j*$Window_size_PAR}{"sperm"} = $accum_value;
}

##draw HapMapII recombination rate curve, the value is cM/Mb
$accum_value = 0.0;
my $Hapmap_p = $HapMapRate{$chr};
for (my $j=1; $j<$Window_num; $j++) {
	$Hapmap_p->[$j-1] ||= 0;
	$Hapmap_p->[$j] ||= 0;
	my $rate1 = ($Hapmap_p->[($j-2>=0) ? ($j-2) : 0] + $Hapmap_p->[$j-1] + $Hapmap_p->[$j] ) / (3*$Window_size_PAR / 1000000) / $Max_cM_vs_Mb_ratio;
	$rate1 = 1 if($rate1 > 1);
	my $x1 = $x_start + $rate1*$x_width;
	my $y1 = $top_edge_pxi + ($j-1)*$win_y_shift;
	my $rate2 = ($Hapmap_p->[$j-1] + $Hapmap_p->[$j] +$Hapmap_p->[($j+1<$Window_num) ? ($j+1) : $j]) / (3*$Window_size_PAR / 1000000) / $Max_cM_vs_Mb_ratio;
	$rate2 = 1 if($rate2 > 1);
	my $x2 = $x_start + $rate2*$x_width;
	my $y2 = $top_edge_pxi + $j*$win_y_shift;
	$svg->line('x1',$x1,'y1',$y1,'x2',$x2,'y2',$y2,'stroke','orange','stroke-width',1);
	
	$accum_value += $Hapmap_p->[$j-1];
	$StatAccum{$chr}{$j*$Window_size_PAR}{"hapmap"} = $accum_value;
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
$svg->text('x',$x1-4,'y',$y2+15,'fill','black', 'font-family','Arial','font-size',12, '-cdata', "0");
$svg->text('x',$x2-9,'y',$y2+15,'fill','black', 'font-family','Arial','font-size',12, '-cdata', "$Max_cM_vs_Mb_ratio");

my $svg_file = "recombination_rate_sperms_HapMap_deCODE.win$Window_size.3Times.svg";
open OUT, ">$svg_file" || die "fail create $svg_file\n";
print OUT $svg->xmlify();
close OUT;
`svg2xxx.pl  -memory 10G  $svg_file`;


#exit;

##draw accumulative cM vs Mb curves
mkdir("chr_accumu_Mb_vs_cM_curves");
for (my $i=1; $i<=23; $i++) {
	my $chr = "chr$i";
	$chr =~ s/23/X/g;
	my $lst_file = "chr_accumu_Mb_vs_cM_curves/$chr.accumulative.Mb.cM.lst";
	my $win_size = ($chr eq "chrX") ? $Window_size_PAR : $Window_size;
	my $Window_num = int($ChrLen{$chr} / $win_size) + 1;
	my $Xaxis_len = ($chr ne "chrX") ? $ChrLen{$chr} : $ChrLen{$chr}*20;
	$Xaxis_len = 50*(int($Xaxis_len / 50000000) + 1);
	my $Yaxis_len = $Xaxis_len * 1.3;
	my $X_unit = ($chr eq "chrX") ? "1/20 Mb" : "Mb";
	my $note = ($chr eq "chrX") ? "par1" : $chr;
	my $head = "Type:Line
LineWidth:3
Width:640
Height:480
WholeScale:0.8
MarkPos:tl
MarkScale:0.8
MarkNoBorder:0
MarkStyle:v
FontSize:46
FontFamily:ArialNarrow-Bold
Note:$note
X:Physical length ($X_unit)
Y:Genetic distance (cM)
Xstart:0
Xstep:50 
Xend:$Xaxis_len
XCut:1
Ystart:0
Ystep:50
Yend:$Yaxis_len
YCut:1
XScaleDiv:10
YScaleDiv:10
Note2:
:End
:End
";
	my ($output1,$output2,$output3);
	$output1 .= "\n\nColor: blue\nMark: Sperm\n";
	$output2 .= "\n\nColor: orange\nMark: HapMapII\n";
	$output3 .= "\n\nColor: cyan\nMark: deCODEmale\n";
	for (my $j=1; $j<$Window_num; $j++) {
		my $accum_len = $j*$win_size;
		my $accum_len_Mb = ($chr ne "chrX") ? int($accum_len / 1000000) : int($accum_len / 50000);
		$output1 .= "$accum_len_Mb: ".$StatAccum{$chr}{$accum_len}{"sperm"}."\n";
		$output2 .= "$accum_len_Mb: ".$StatAccum{$chr}{$accum_len}{"hapmap"}."\n";
		$output3 .= "$accum_len_Mb: ".$StatAccum{$chr}{$accum_len}{"decode"}."\n";
		
	}

	open OUT, ">$lst_file" || die "fail";
	print OUT $head.$output1.$output2.$output3;
	close OUT;
	`/gpfsdata/Analysis/fanwei/bin/distribute_svg/distribute_svg.pl  $lst_file  $lst_file.svg`;
	`svg2xxx.pl $lst_file.svg`;

}

##calculate the correlation between sperm&hapmap, sperm&decode, decode&hapmap;
open OUT, ">correlation_between_genetic_maps.tab" || die "fail";
print OUT "#Chr\tSperm&HapMap\tSperm&deCODE\tHapMap&deCODE\tWinNum\n";
my (@spermmap, @hapmap, @decodemap);
for (my $i=1; $i<=22; $i++) {
	my $chr = "chr$i";
	$chr =~ s/23/X/g;
	my $sperm_p = $StatCorre{$chr}{"sperm"};
	my $hapmap_p = $StatCorre{$chr}{"hapmap"};
	my $decode_p = $StatCorre{$chr}{"decode"};
	my $wins_num = @$sperm_p;
	push @spermmap, @$sperm_p;
	push @hapmap, @$hapmap_p;
	push @decodemap, @$decode_p;
	my $r_sperm_hapmap = regression_LSM($sperm_p, $hapmap_p);
	my $r_sperm_decode = 0;
	my $r_hapmap_decode = 0;
	if ($chr ne "chrX") {
		$r_sperm_decode = regression_LSM($sperm_p, $decode_p);
		$r_hapmap_decode = regression_LSM($hapmap_p, $decode_p);
	}
	print OUT "$chr\t$r_sperm_hapmap\t$r_sperm_decode\t$r_hapmap_decode\t$wins_num\n";
}
my $wins_num = @spermmap;
my $r_sperm_hapmap = regression_LSM(\@spermmap, \@hapmap);
my $r_sperm_decode = regression_LSM(\@spermmap, \@decodemap);
my $r_hapmap_decode = regression_LSM(\@hapmap, \@decodemap);
print OUT "Total\t$r_sperm_hapmap\t$r_sperm_decode\t$r_hapmap_decode\t$wins_num\n";
close OUT;



####################################################
################### Sub Routines ###################
####################################################


##caculate regression equation by the least square method: y = ax + b
sub regression_LSM()
{	
	my ($p1, $p2) = @_;
	my @x = @$p1;
	my @y = @$p2;
	my ($a, $b, $r);
	my $N = @x;
	my $Mean_X=0.0;
	my $Mean_Y=0.0;
	my $S_XX=0.0;
	my $S_XY=0.0;
	my $S_YY=0.0;

	for (my $i=0; $i<$N; $i++){
		$Mean_X += $x[$i];
		$Mean_Y += $y[$i];
	}
	$Mean_X /= $N;
	$Mean_Y /= $N;

	for (my $i=0; $i<$N; $i++){
		$S_XX +=  ($x[$i] - $Mean_X) * ($x[$i] - $Mean_X);
		$S_XY +=  ($x[$i] - $Mean_X) * ($y[$i] - $Mean_Y);
		$S_YY +=  ($y[$i] - $Mean_Y) * ($y[$i] - $Mean_Y);
	}

	$a = $S_XY / $S_XX;
	$b = $Mean_Y - $a * $Mean_X;
	$r = $S_XY / sqrt($S_XX * $S_YY);
	return $r;
}


sub read_chrlen{
	my $file = shift;
	my $ChrLen_p = shift;

	open IN, $file || die "fail $file";
	while (<IN>) {
		next if(/^\#/);
		chomp;
		my ($chr,$len) = ($1,$2) if(/^(\S+)\s+(\d+)/);
		$len = 2639520 if($chr eq "chrX");
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
		my $win_size = ($chr ne "chrX") ? $Window_size : $Window_size_PAR;
		#print "$chr\t$Window_size\t$Window_size_PAR\t$win_size\n";
		if ($chr eq "chrX") {
			if ($pos > 2699520 ) {
				next;
			}
		}
		my $bases_num = $t[4];
		my $win_id = int($pos / $win_size);
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
		#$g->line('x1',$X,'y1',$rulcfg{Y} - $scale_size,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '1');
		#$g->line('x1',$X,'y1',$rulcfg{Y} + $scale_size,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1)  if ($rulcfg{rulerstyle} eq '2');
		#$g->line('x1',$X,'y1',$rulcfg{Y} - $scale_size,'x2',$X,'y2',$rulcfg{Y},'stroke','#000000','stroke-width',1) if ($rulcfg{rulerstyle} eq '3');
		my $disp_scale_text = $i / $rulcfg{bigscalesize} if(defined $rulcfg{bigscalesize});
		#$rulcfg{svg}->text('x',$X - textWidth($rulcfg{font_family},$rulcfg{font_size},$disp_scale_text) / 2,'y',$rulcfg{Y}+textHeight($rulcfg{font_size}),'fill','#000000','-cdata',$disp_scale_text,'font-size',$rulcfg{font_size}, 'font-family',$rulcfg{font_family}) if ($rulcfg{rulerstyle} eq '1');
		#$rulcfg{svg}->text('x',$X - textWidth($rulcfg{font_family},$rulcfg{font_size},$disp_scale_text) / 2,'y',$rulcfg{Y}+$scale_size+textHeight($rulcfg{font_size})+2,'fill','#000000','-cdata',$disp_scale_text,'font-size',$rulcfg{font_size}, 'font-family',$rulcfg{font_family})  if ($rulcfg{rulerstyle} eq '2');
		#$rulcfg{svg}->text('x',$X - textWidth($rulcfg{font_family},$rulcfg{font_size},$disp_scale_text) / 2,'y',$rulcfg{Y}-$scale_size-textHeight($rulcfg{font_size})+6,'fill','#000000','-cdata',$disp_scale_text,'font-size',$rulcfg{font_size}, 'font-family',$rulcfg{font_family}) if ($rulcfg{rulerstyle} eq '3');
	}

}

##read not used sperm list
sub read_sperm_list{
	my $file = shift;
	open IN, $file || die "$file";
	while (<IN>) {
		next if(/^\#/);
		my $sperm = $1 if(/(S\d\d)/);
		$NotUsedSperms{$sperm} = 1;
	}
	close IN;
}


sub read_cross_over{
	my $file = shift;

	open IN, $file || die "fail ";
	while (<IN>) {
		next if(/^\#/);
		chomp;
		my @t = split /\t/;
		my $sperm_id = $t[0];
		my $chr = $t[1];
		my $pos = $t[2];
		my $win_size = ($chr ne "chrX") ? $Window_size : $Window_size_PAR;
		if ($chr eq "chrX") {
			#$pos -= 150000000 if($pos >= 3000000);
			next if($pos > 2639520);
		}
		my $status = $t[6];
		$status =~ s/\s//g;

		if(exists $NotUsedSperms{$sperm_id}){
			print STDERR "skipped $sperm_id\n";
			next;
		}
		
		if ( $status == 1 ) {
			push @{$CrossOver{$chr}}, $pos;
			$StatRecomb{$sperm_id}{$chr} ++;
			my $win_id = int($pos / $win_size);
			$SpRecombRate{$chr}[$win_id] ++;
		} 
		
	}
	close IN;

}


sub read_HapMap_file {
	my $file = shift;
	
	my %data;
	open IN, $file || die "fail ";
	while (<IN>) {
		next if(/^\#/);
		chomp;
		my @t = split /\s+/;
		my $chr = $t[0];
		my $pos = $t[1];
		
		##&& $pos > 1029793 && $pos < 1459643
		if ($chr eq "chrX" ) {
			print "#$chr\t$pos\t$t[2]\t$t[3]\n";
		}

		if ($chr eq "chrX") {
			#$pos -= 150000000 if($pos >= 3000000);
			next if($pos > 2639520);
		}
		my $cM = $t[3];
		$data{$chr}{$pos} = $cM;
	}
	close IN;

	foreach my $chr (sort keys %data) {
		my $chr_p = $data{$chr};
		my $last_cM = 0.0;
		my $win_size = ($chr ne "chrX") ? $Window_size : $Window_size_PAR;
		foreach my $pos (sort {$a<=>$b} keys %$chr_p) {
			my $cM = $chr_p->{$pos} - $last_cM;
			if ($chr_p->{$pos} < $last_cM) {
				print STDERR "Error happens $chr_p->{$pos}  $last_cM\n";
			}
			$last_cM = $chr_p->{$pos};
			my $win_id = int($pos / $win_size);
			if ($chr ne "chrX") {
				$HapMapRate{$chr}[$win_id] += $cM;
			}else{ ##because HapMap is averaged-sex cM, but PAR only exists in male
				$HapMapRate{$chr}[$win_id] += 2*$cM;
			}
			
			$StatHapMap{$chr} = $last_cM;
		}
	}
	$StatHapMap{"chrX"} *= 2;
	my $total_cM = 0.0;
	foreach my $chr (sort keys %StatHapMap) {
		$total_cM += $StatHapMap{$chr};
	}
	$StatHapMap{"Total"} = $total_cM;
}


##chr1    rs1181875       3681831 0.053508
sub read_deCODE_file{

	my $file = shift;
	
	my %data;
	open IN, $file || die "fail ";
	while (<IN>) {
		next if(/^\#/);
		chomp;
		my @t = split /\t/;
		my $chr = $t[0];
		my $win_size = ($chr ne "chrX") ? $Window_size : $Window_size_PAR;
		my $pos = $t[2];
		if ($chr eq "chrX") {
			#$pos -= 150000000 if($pos >= 3000000);
			next if($pos > 2639520);
		}
		my $cM = $t[3];
		my $win_id = int($pos / $win_size);
		$deCodeRate{$chr}[$win_id] += $cM;
		$StatDeCode{$chr} += $cM;
	}
	close IN;

	my $total_cM = 0.0;
	foreach my $chr (sort keys %StatDeCode) {
		$total_cM += $StatDeCode{$chr};
	}
	$StatDeCode{"Total"} = $total_cM;

	##Note that there is no PAR data in deCODE data.

}

