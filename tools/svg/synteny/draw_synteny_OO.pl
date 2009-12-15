#!/bin/env perl

=head1 Name

draw_synteny_OO.pl  --  draw synteny figure between two sequences (OO version)

=head1 Description

draw linear figures for synteny analysis with nucmer result,
also add gene and TE elements to svg figure.

=head1 Version

  Author: Hu Xuesong, huxuesong@genomics.org.cn
  Version: 1.0,  Date: 2009-12-14

=head1 Usage

  perl draw_synteny_parallel_melon4.pl [options]
  --width <int>         set the image width (1200)
  --height <int>        set the image height (720)
  --chr <file>          set the chromosome length file
  --list <file>         set the pairwise relation list file
//  --scale <int>         set the unit for big scales
//  --resolution <float>  set the figure resolution, default 100000
  --synteny_heigth <int>  set the height of each unit of synteny figure, default=60
  --gene1 <file>         draw exon-intron structure for one gene set
  --gene2 <file>         draw exon-intron structure for another gene set
  --TE1 <file>           draw colorful blocks for one TE set
  --TE2 <file>           draw colorful blocks for another TE set
  --contig1 <file>		 draw contigs for one scaffold set
  --contig2	<file>		 draw contigs for another scaffold set
  --outdir <str>        set the output directory
  --verbose   output verbose information to screen
  --help      output help information to screen

=head1 Exmple

  perl draw_synteny_OO.pl  --synteny_heigth 150 --scale 10000 --resolution 100 -chr ./melon_cucumber_chrs.list --list ./melo_4BACs.delta.filter.coords.list --gene1 ./melo_4BACs.gbk.noTE.gff -TE1 ./melo_4BACs.fa.TE.gff --gene2 ./cucumber.3scaffs.gene.gff.need.noTE.gff --TE2 ./cucumber.scaffs.TE.gff

=cut

use strict;
use warnings;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use Data::Dumper;
use lib "$Bin/../lib";
use lib 'U:\dev\sda\5\develop\ActivePerl510\Perl64\site\lib\\';
use SVG;
use FontSize;

#my $figure_resolution; ##one point stands for 10000 bp
#my $big_scale;
my ($synteny_heigth);
my ($Gene1_file,$TE1_file,$Gene2_file,$TE2_file,$contig1_file,$contig2_file,$depth_file,$img_width,$img_height);
my ($pairwise_file,$chr_len_file1);
my ($Verbose,$Help,$Outdir);
GetOptions(
	"width:i"=>\$img_width,
	"height:i"=>\$img_height,
	"chr:s"=>\$chr_len_file1,
	"list:s"=>\$pairwise_file,
#	"scale:i"=>\$big_scale,
#	"resolution:f"=>\$figure_resolution,
	"synteny_heigth:i"=>\$synteny_heigth,
	"gene1:s"=>\$Gene1_file,
	"gene2:s"=>\$Gene2_file,
	"TE1:s"=>\$TE1_file,
	"TE2:s"=>\$TE2_file,
	"contig1:s"=>\$contig1_file,
	"contig2:s"=>\$contig2_file,
	"outdir:s"=>\$Outdir,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
$img_width ||=1200;
$img_height ||=720;
#$figure_resolution ||= 100000;
#$big_scale ||= 10000000;
$synteny_heigth ||= 60;
$Outdir ||= ".";
die `pod2text $0` if (!$chr_len_file1 || !$pairwise_file || $Help);

#Contig_4	1	97799	scaffold375	1	34612	-
#Contig_4	1	97799	scaffold2282	377	67038	+
sub read_chromosome_file($$) {
	my $chromosome_file=shift;
	my $chromosome_p=shift;
	open CHR,$chromosome_file or die "[x]$!";
	while(<CHR>) {
		chomp;
		next if $_=~/^chr_id\t/;
		s/^\s+//;
		next if(!$_);
		my @t = split /\s+/;
		push @{$chromosome_p->{$t[0]}}, \@t;
	}
	close CHR;
}

#Contig_1_4	Contig_1	-	79705	59416	scaffold629_4	scaffold629	+	2529462	2550090
#Contig_2_5	Contig_2	+	1	6619	scaffold421_5	scaffold421	+	290060	296703
sub read_pairwise_file {
	my $pairwise_file = shift;
	my $synteny_p = shift;
	open IN, $pairwise_file || die "[x]Fail open $pairwise_file !";
	while (<IN>) {
		chomp;
		my ($gene1,$chr1,$strand1,$start1,$end1,$gene2,$chr2,$strand2,$start2,$end2) = split /\s+/;
		push @$synteny_p,[$gene1,$start1,$end1,$strand1,$gene2,$start2,$end2,$strand2,$chr1,$chr2];
	}
	close IN;
}

#Contig_8	BGI	mRNA	10453	11903	.	+	.	ID=Apl_00006;
#Contig_8	BGI	CDS	10453	10464	.	+	0	Parent=Apl_00006;
#Contig_8	BGI	CDS	11887	11903	.	+	0	Parent=Apl_00006;
sub read_gene_gff($$) {
	my $file=shift;
	my $ref=shift;
	open (IN,$file) || die ("[x]Fail open $file !\n");
	while (<IN>) {
		chomp;
		s/^\s+//;
		next if(/^\#/);
		my @t = split(/\t/);
		my $tname = $t[0];
		my $qname;
		if ($t[2] eq 'mRNA' || $t[2] eq 'CDS'|| $t[2] eq 'CDs') {
			$qname = $1 if($t[8] =~ /^GenePrediction\s+(\S+)/ || $t[8] =~ /^ID=(\S+?);/ || $t[8] =~ /^Parent=(\S+?);/);
		}
		if ($t[2] eq 'match' || $t[2] eq 'HSP') {
			$qname = $1 if($t[8] =~ /Target\s+\"(\S+)\"/);
		}
		if ($t[2] eq 'mRNA' || $t[2] eq 'match') {
			my $start = $t[3];
			my $end = $t[4];
			my $strand = $t[6];
			$ref->{$tname}{$qname}{strand} = $strand;
			$ref->{$tname}{$qname}{start} = $start;
			$ref->{$tname}{$qname}{end} = $end;
		}
		if ($t[2] eq 'CDS' || $t[2] eq 'HSP'|| $t[2] eq 'CDs') {
#			print "$tname\t$qname\n";
			push @{$ref->{$tname}{$qname}{exon}}, [$t[3],$t[4]] if(exists $ref->{$tname}{$qname});
		}
	}
	close(IN);

	foreach my $scaff (sort keys %$ref) {
		my $scaff_p = $ref->{$scaff};
	#	print Dumper ($ref->{$scaff});
		foreach my $gene (sort keys %$scaff_p) {
		#	print Dumper ($ref->{$scaff});
			my $gene_p = $scaff_p->{$gene};
			my @exon = @{$gene_p->{exon}};
			@exon = reverse @exon if($exon[0][0] > $exon[-1][0]);
			my @new;
			foreach my $p (@exon) {
				push @new,$p->[0],$p->[1];
			}
			$ref->{$scaff}{$gene}{exon} = \@new;
		}
	}
}

##gff-version 3
#scaffold375	RepeatMasker	Transposon	1	3000	1611	+	.	ID=TE055206;Target=CR1-Y2_Aves 340 1177;Class=LINE/CR1;PercDiv=29.0;PercDel=2.9;PercIns=0.2;
sub read_TE_gff($$) {
	my $file=shift;
	my $ref=shift;
	open (IN,$file) || die ("[x]Fail open $file !\n");
	while (<IN>) {
		chomp;
		s/^\s+//;
		next if(/^\#/);
		my @t = split(/\t/);
		my $tname = $t[0];
		my $qname = $1 if($t[8] =~ /^ID=([^;]+);*/);
		my $type = $1 if($t[8] =~ /Class=([^\/;]+)/);
		my $start = $t[3];
		my $end = $t[4];
		my $strand = $t[6];
		$ref->{$tname}{$qname} = [$start,$end,$strand,$type];
	}
	close(IN);
}

#Contig_4	1	7836
#Contig_4	8956	26732
sub read_contig_list($$) {
	my $file=shift;
	my $hash_p=shift;
	open IN,$file or die "[x]Fail open $file !\n";
	while(<IN>) {
		chomp;
		my @ary=split /\t/;
		push @{$hash_p->{$ary[0]}},[$ary[1],$ary[2]];
	}
	close IN;
}






my $pairwise_file_base = basename($pairwise_file);
$Outdir =~ s/\/$//;
mkdir $Outdir unless (-d $Outdir);

my %Chr1;
my @Synteny;

#$figure_resolution = 1 / $figure_resolution;

my (%SVG,@Show);
my @Order=('Ruler1','mRNA1','TE1','Name1','Bar1','Synteny',
	'Bar2','Name2','TE2','mRNA2','Ruler2');
our %Heights=(
	Ruler1 => 20,
	Ruler2 => 20,
	mRNA1 => 15,
	mRNA2 => 15,
	TE1 => 15,
	TE2 => 15,
	Name1 => 10,
	Name2 => 10,
	Bar1 => 15,
	Bar2 => 15,
	Synteny => $synteny_heigth
);
my $Y_junt=10;
my %toShow=(
	Ruler1 => 1,
	Ruler2 => 1,
	Name1 => 1,
	Name2 => 1,
	Bar1 => 1,
	Bar2 => 1,
	Synteny => 1,
);

##read data into memory
read_chromosome_file($chr_len_file1,\%Chr1);
read_pairwise_file($pairwise_file,\@Synteny);

my (%Gene1,%Gene2,%TE1,%TE2,%contig1,%contig2,%depth);
read_gene_gff($Gene1_file,\%Gene1) if(defined $Gene1_file);
read_gene_gff($Gene2_file,\%Gene2) if(defined $Gene2_file);
read_TE_gff($TE1_file,\%TE1) if(defined $TE1_file);
read_TE_gff($TE2_file,\%TE2) if(defined $TE2_file);
read_contig_list($contig1_file,\%contig1) if(defined $contig1_file);
read_contig_list($contig2_file,\%contig2) if(defined $contig2_file);
##figure parameters setting
our $font = FontSize->new();
our $font_family = "Arial";
#my $font_size = 24;	# stringHeight = int($size*0.72); => font_size = stringHeight / 0.72
our $x_margin = 20;
my $y_margin = 20;
my $Seq2_junt=1000;

my %Synteny_length_bp;
for my $key (keys %Chr1) {
	my ($Seq1_name,$Seq1_start,$Seq1_end) = @{${$Chr1{$key}}[0]};
	my $Seq2_lenA = 0;
	my $Seq1_len = $Seq1_end - $Seq1_start +1;
	for (@{$Chr1{$key}}) {
		my ($Seq1_name,$Seq1_start,$Seq1_end,$Seq2_name,$Seq2_start,$Seq2_end,$Reverse) = @$_;
		#$Seq1_lenA = $Seq1_end - $Seq1_start;
		$Seq2_lenA += $Seq2_junt if $Seq2_lenA;
		$Seq2_lenA += $Seq2_end - $Seq2_start +1;
	}
	$Synteny_length_bp{$key} = ($Seq1_len >= $Seq2_lenA) ? $Seq1_len : $Seq2_lenA;
	my $X_start=0;

	if(defined $contig1_file && exists $contig1{$Seq1_name}) {
		$toShow{'Bar1'}=1;
		for (@{$contig1{$Seq1_name}}) {
			my ($start,$end)=@$_;
			my $str='1: '.$Seq1_name;
			push @{$SVG{$key}{'Bar1'}},[$str,$start-$Seq1_start+1,$end-$start+1,'black'];#($str,$margin+$start,$len,$color)
		}
	}

	$SVG{$key}{'Name1'}=[[$Seq1_name,0,$Seq1_len]];
	if (defined $TE1_file && exists $TE1{$Seq1_name}) {
		$toShow{'TE1'}=1;
		my $Seq1_p = $TE1{$Seq1_name};
		foreach my $TE1_id (sort keys %$Seq1_p) {
			my $TE1_p = $Seq1_p->{$TE1_id}; ##$ref->{$tname}{$qname} = [$start,$end,$strand,$type];
			my ($start,$end,$strand,$type) = @$TE1_p;
			my $str=join(': ',$type,$TE1_id);
			push @{$SVG{$key}{'TE1'}},[$str,$start-$Seq1_start+1,$end-$start+1];#($str,$margin+$start,$len)
			}
	}

	if (defined $Gene1_file && exists $Gene1{$Seq1_name}) {
		$toShow{'mRNA1'}=1;
		my $Seq1_p = $Gene1{$Seq1_name};
			foreach my $gene1_id (keys %$Seq1_p) {
				my @StartEndStrandFull=($Seq1_p->{$gene1_id}{start}-$Seq1_start+1,$Seq1_p->{$gene1_id}{end}-$Seq1_start+1,$Seq1_p->{$gene1_id}{strand},1);
				#next if ($StartEndStrand[0] > $Seq1_end or $StartEndStrand[1] < $Seq1_start);
				if ($StartEndStrandFull[0] < $Seq1_start-$Seq2_junt/2) {
					$StartEndStrandFull[3]=0;
					$StartEndStrandFull[0]=$Seq1_start-$Seq2_junt/2;
				}
				if ($StartEndStrandFull[1] > $Seq1_end+$Seq2_junt/2) {
					$StartEndStrandFull[3]=0;
					$StartEndStrandFull[1]=$Seq1_end+$Seq2_junt/2;
				}
				my @exon = @{$Seq1_p->{$gene1_id}{exon}};
				#my $gene1_strand = $Seq1_p->{$gene1_id}{strand};
				my @X_exon;
				#push @X_exon,$_-$Seq1_start+1 for (@exon);
				for (@exon) {
					my $pos=$_-$Seq1_start+1;
					$pos=-$Seq2_junt/2 if $pos<-$Seq2_junt/2;
					$pos=$Seq1_end+$Seq2_junt/2 if $pos>$Seq1_end+$Seq2_junt/2;
					push @X_exon,$pos;
				}
				push @{$SVG{$key}{'mRNA1'}},[$gene1_id,\@StartEndStrandFull,\@X_exon];
			}
		}

	for (@{$Chr1{$key}}) {
		my ($Seq1_name,$Seq1_start,$Seq1_end,$Seq2_name,$Seq2_start,$Seq2_end,$Reverse) = @$_;
		#my $Seq1_len = $Seq1_end - $Seq1_start +1;
		my $Seq2_len = $Seq2_end - $Seq2_start +1;
		#$SVG{$key}{'Name1'}=[[$Seq1_name,0,$Seq1_len]];
		push @{$SVG{$key}{'Name2'}},[$Seq2_name,$X_start,$Seq2_len];
		if (defined $TE2_file && exists $TE2{$Seq2_name}) {
			$toShow{'TE2'}=1;
			my $Seq2_p = $TE2{$Seq2_name};
			foreach my $TE2_id (sort keys %$Seq2_p) {
				my $TE2_p = $Seq2_p->{$TE2_id}; ##$ref->{$tname}{$qname} = [$start,$end,$strand,$type];
				my ($start,$end,$strand,$type) = @$TE2_p;
				next if ($end<$Seq2_start or $start>$Seq2_end);
				my $str=join(': ',$type,$TE2_id);
				my $len=$end-$start+1;
				if ($Reverse eq "-") {
					$start = $Seq2_end - $end +1;
				}
				my $begin=$start-$Seq2_start+1;
				if ($begin<-$Seq2_junt/3) {
					$len += $begin+$Seq2_junt/3;
					$begin=-$Seq2_junt/3;
				}
				push @{$SVG{$key}{'TE2'}},[$str,$X_start+$begin,$len];
			}

		}


		$X_start += $Seq2_len + $Seq2_junt;
	}
}

our $WidthPoint;
our $HeightPoint;
for (@Order) {
	if (defined $toShow{$_}) {
		push @Show,$_;
		$HeightPoint += $Heights{$_} + $Y_junt;
	}
}
$HeightPoint -= $Y_junt;
$HeightPoint = ($img_height-2*$y_margin)/$HeightPoint;
our ($Y_start,$svg);

for my $key (keys %SVG) {
	$WidthPoint = ($img_width-2*$x_margin)/$Synteny_length_bp{$key};
	$svg = SVG->new('width',$img_width,'height',$img_height);
$svg->rect('x',0,'y',$img_height/2,'height',$y_margin,'width',$img_width,'fill','black');
$svg->rect('x',$x_margin,'y',$y_margin,'height',$img_height-2*$y_margin,'width',$img_width-2*$x_margin,'stroke','black','fill','rgb(196,196,196)');
	$Y_start=$y_margin;
	open OUT, ">$Outdir/$key.synteny.svg" || die "[x]Fail create $!";
	for my $ss (@Show) {
		no strict 'refs';
		$ss =~ /^([^\d]+)(\d)$/;
		my $Func='SVG_'.$1;
		my $OE=($2==2);
		&$Func($ss,$_,$OE) for (@{$SVG{$key}{$ss}});
		$Y_start += ($Y_junt + $Heights{$ss})*$HeightPoint;
	}
	$svg->rect('x',0,'y',0,'height',$img_height,'width',$x_margin/2,'fill','white');
	$svg->rect('x',$img_width-$x_margin/2,'y',0,'height',$img_height,'width',$x_margin/2,'fill','white');
	print OUT $svg->xmlify();
	close OUT;
}
print "$WidthPoint,$HeightPoint,$svg\n";
sub SVG_Name($$$) {
$svg->rect('x',$x_margin, 'y',$Y_start,'width',$img_width-2*$x_margin,'height',$Heights{$_[0]}*$HeightPoint,'fill','none','stroke','blue');
	my ($str,$margin,$len)=@{$_[1]};
	my $height=$Heights{$_[0]}*$HeightPoint;#*.8;
	#my $OE=$_[2];
	my $Y0=$Y_start+$_[2]*$Heights{$_[0]}*$HeightPoint;#*.8;
	my $font_size = int($height / .72);
	my $name_width = $font->stringWidth($font_family,$font_size,$str);
#$name_width=0;
	my $x=$x_margin+(($len+$margin+$margin)*$WidthPoint-$name_width)/2;
	$svg->text('x',$x,'y',$Y0+$height*(!$_[2]),'-cdata',$str,'font-family',$font_family,'font-size',$font_size);
#$svg->line('x1',$x_margin,'y1',$Y_start,'x2',$x+$name_width/2,'y2',$Y0,'stroke','#000000','stroke-width',1);
	print "$_[0],$str,$margin,$len\t$Y_start\t$name_width,$font_size\n";
}

sub SVG_Ruler($$) {
$svg->rect('x',$x_margin, 'y',$Y_start,'width',$img_width-2*$x_margin,'height',$Heights{$_[0]}*$HeightPoint,'fill','none','stroke','blue');
}
sub SVG_mRNA($$) {
$svg->rect('x',$x_margin, 'y',$Y_start,'width',$img_width-2*$x_margin,'height',$Heights{$_[0]}*$HeightPoint,'fill','none','stroke','blue');
#[$gene1_id,\@StartEndStrand,\@X_exon];
	my ($str,$SESFrefA,$EXONrefA)=@{$_[1]};
	my $height=$Heights{$_[0]}*$HeightPoint;
	my $M_Height=$Heights{$_[0]}*.824;
	$M_Height = $Y_junt if $M_Height > $Y_junt;
	$M_Height *= .75;
	my $font_size = int($M_Height*$HeightPoint / .72);
	my $M_str=($$SESFrefA[3] == 1)?'M':'m';
	my $M_width = $font->stringWidth($font_family,$font_size,$M_str);
	my $Y0=$Y_start+( (!$_[2])*$Heights{$_[0]} + $M_Height*2*(.5-$_[2]) )*$HeightPoint;
	my $g = $svg->group('id'=>$str,'onclick',"alert('$str')");
	# Bone Line
	my $y=$Y_start+$height/2;
	my $x1=$x_margin+$$SESFrefA[0]*$WidthPoint;
	my $x2=$x_margin+$$SESFrefA[1]*$WidthPoint;
	$g->line('x1',$x1,'y1',$y,'x2',$x2,'y2',$y,'stroke','rgb(0,96,196)','stroke-width',$height/4);
	# Exon Boxes
	while (@$EXONrefA > 1) {
		my $x1=shift @$EXONrefA;
		my $x2=shift @$EXONrefA;
		my $w=$x2-$x1;
		$x1 *= $WidthPoint;	$w *= $WidthPoint;
		$w=1 if $w<1;
		$x1 += $x_margin;
		$g->rect('x',$x1, 'y',$Y_start,'width',$w,'height',$height,'fill','rgb(0,96,96)');
	}
	# M
	my $A_PointX;
	if ($$SESFrefA[2] eq '+') {
		$A_PointX=$x2;
		$A_PointX=$img_width-$x_margin/2 if $A_PointX>$img_width-$x_margin/2;
		$A_PointX -= $M_width;
	} else {
		$A_PointX=$x1;
		$A_PointX=$x_margin/2 if $A_PointX<$x_margin/2;
	}
	$g->text('x',$A_PointX,'y',$Y0,'-cdata',$M_str,'font-family',$font_family,'font-size',$font_size);
}
sub SVG_TE($$) {
$svg->rect('x',$x_margin, 'y',$Y_start,'width',$img_width-2*$x_margin,'height',$Heights{$_[0]}*$HeightPoint,'fill','none','stroke','blue');
	my ($str,$start_margin,$len)=@{$_[1]};
	my $height=$Heights{$_[0]}*$HeightPoint;
	my $x=$x_margin+$start_margin*$WidthPoint;
	my $w=$len*$WidthPoint;
	$w=1 if $w<1;
	$svg->rect('x',$x, 'y',$Y_start,'width',$w,'height',$height,'fill','red','onclick',"alert('$str')");
}
sub SVG_Bar($$) {
$svg->rect('x',$x_margin, 'y',$Y_start,'width',$img_width-2*$x_margin,'height',$Heights{$_[0]}*$HeightPoint,'fill','none','stroke','blue');
#($str,$margin+$start,$len,$color)
	my ($str,$start_margin,$len,$color)=@{$_[1]};
	my $height=$Heights{$_[0]}*$HeightPoint;
	my $x=$x_margin+$start_margin*$WidthPoint;
	my $w=$len*$WidthPoint;
	$w=1 if $w<1;
	$svg->rect('x',$x, 'y',$Y_start,'width',$w,'height',$height,'fill',$color,'onclick',"alert('$str')");
}
sub SVG_Synteny($$) {
$svg->rect('x',$x_margin, 'y',$Y_start,'width',$img_width-2*$x_margin,'height',$Heights{$_[0]}*$HeightPoint,'fill','none','stroke','blue');
}
