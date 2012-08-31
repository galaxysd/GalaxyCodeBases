#!/usr/bin/perl

=head1 Name



=head1 Description



=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6
  Note:

=head1 Usage
  perl draw_sperm.pl  <chr_len_file> <chr_gap_file> <chr_centromere_file> <sperms_halotype_dir> 
  --window <int>   window size, i.e number of genomic bases in each window, default=1000000
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
my $cross_over_file = shift;

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
read_cross_over($cross_over_file);

my $left_edge_pxi = 50;
my $right_edge_pxi = 20;
my $top_edge_pxi = 40;
my $bottom_edge_pxi = 30;


my $chr_x_shift = 12;
my $win_y_shift = 1;
my $chr_width = $chr_x_shift*0.6;

for (my $c=1; $c<=$Chrom_num; $c++) {
	my $chr = "chr$c";
	$chr =~ s/23/X/;
	$chr =~ s/24/Y/;
	next if($chr eq "chrX" || $chr eq "chrY");
	
	my $Max_win_num = int($ChrLen{$chr} / $Window_size) + 1; ##chr1 is the largest chromosome
	my $figure_width = $left_edge_pxi + $right_edge_pxi + $Sperm_num*$chr_x_shift;
	my $figure_height = $top_edge_pxi + $bottom_edge_pxi + $win_y_shift*$Max_win_num;

	my $svg = SVG->new('width',$figure_width,'height',$figure_height);
	$svg->rect('x',0, 'y',0,'width',$figure_width,'height',$figure_height,'fill',"white");

	$svg->text('x',$left_edge_pxi-40,'y',$figure_height/2,'fill','black', 'font-family','Arial','font-size',16, '-cdata',$chr);

	for (my $i=0; $i<$Sperm_num; $i++) {
		my $k = $i + 1;
		my $sperm_id = ($i < 9) ? "S0$k" : "S$k";

		my $x = $left_edge_pxi + $i*$chr_x_shift;
		my $y;
		my $text_x = $x + 8;
		my $text_y = $top_edge_pxi - 5;
		my $transform_format = "rotate(-90 $text_x,$text_y)";
		$svg->text('x',$text_x ,'y',$text_y,'fill','black', 'font-family','Arial','font-size',12, '-cdata',$sperm_id, 'transform',$transform_format);
		my $Window_num = int($ChrLen{$chr} / $Window_size) + 1;
		
		##draw the automatic identified crossover sites, as background
		my $cross_p = $CrossOver{$sperm_id}{$chr};
		my $status_p = $CrossOverStatus{$sperm_id}{$chr};
		foreach my $pos (@$cross_p) {
			my $cy = $top_edge_pxi + int($pos/$Window_size)*$win_y_shift;
			my $cx = $x+$chr_width/2;
			my $r = $chr_width*1.1;
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
		my $x = $left_edge_pxi + $i*$chr_x_shift;
		my $Sp_Chr_p = $SpColor{$chr}[$i];
		for (my $win_id = 0; $win_id < $Window_num; $win_id++) {
			my $p = $Sp_Chr_p->[$win_id];
			$p->[0] ||= 0;  ##if it is empty, set default value 0
			$p->[1] ||= 0;  ##if it is empty, set default value 0
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
			
			$y = $top_edge_pxi + $win_id*$win_y_shift;
			$svg->rect('x',$x, 'y',$y,'width',$chr_width,'height',$win_y_shift,'stroke',"none", 'fill',$color);
			
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
		}

		##draw chromosome frame
		$y = $top_edge_pxi;
		$svg->rect('x',$x, 'y',$y,'width',$chr_width,'height',$Window_num*$win_y_shift,'stroke',"black", "fill", "none");


		##draw centromere (draw on the top of figure)
		my ($centro_start,$centro_end) = ($1,$2) if($Centro{$chr} =~ /(\d+),(\d+)/);
		my $centro_pos = int(($centro_start + $centro_end) / 2);
		$y = $top_edge_pxi + int($centro_pos/$Window_size)*$win_y_shift;
		$svg->circle('cx',$x+$chr_width/2,'cy',$y,'r',$chr_width*0.6,'stroke','black','fill','none','stroke-width',3);

	}

	my $svg_file = "$chr.allsperms.haplotype.svg";
	open OUT, ">$svg_file" || die "fail create $svg_file\n";
	print OUT $svg->xmlify();
	close OUT;
	`svg2xxx.pl  -memory 10G  $svg_file`;

}

mkdir("chr-allsperms-haplotype-svg-win$Window_size");
mkdir("chr-allsperms-haplotype-png-win$Window_size");
`mv chr*.allsperms.haplotype.svg  chr-allsperms-haplotype-svg-win$Window_size`;
`mv chr*.allsperms.haplotype.png  chr-allsperms-haplotype-png-win$Window_size`;

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


sub read_cross_over(){
	my $file = shift;

	open IN, $file || die "fail ";
	while (<IN>) {
		next if(/^\#/);
		chomp;
		my @t = split /\s+/;
		my $sperm_id = $t[0];
		my $chr = $t[1];
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

