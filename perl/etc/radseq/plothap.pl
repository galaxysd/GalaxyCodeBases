#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Data::Dump qw(ddx);
use Galaxy;

die "Usage: $0 <marker dat> <tassel dump> [out midfix]\n" if @ARGV<2;
my $markf=shift;
my $inf=shift;
my $outf=shift;

if ($outf) {
	$outf = ".$outf";
} else {
	$outf = '';
}

my $scaffnfo = '/bak/seqdata/2012/tiger/120512_TigerRefGenome/chr.nfo';
my $scaf2locuslst = '/share/users/huxs/work/tiger/20120910/tighap.lst';
print "From [$markf][$inf] to [$inf$outf.(dat|svg)]\n";

my %Stat;

my @Scaffolds = ([qw(scaffold75 scaffold1458)], ['scaffold188']);
my $t=0;
my %scaffolds = map {map {$_ => ++$t} @$_} @Scaffolds;
my %scaff2chr = map { $scaffolds{$_} => $_ } keys %scaffolds;
ddx \@Scaffolds;
ddx \%scaffolds;
ddx \%scaff2chr;
# plothap.pl:32: [["scaffold75", "scaffold1458"], ["scaffold188"]]
# plothap.pl:33: { scaffold1458 => 2, scaffold188 => 3, scaffold75 => 1 }
# plothap.pl:34: { 1 => "scaffold75", 2 => "scaffold1458", 3 => "scaffold188" }

my ($maxLen,@LineLen,%ScaffoldLen,%ScaffoldsOffect)=(0);
open I,'<',$scaffnfo or die $!;
while (<I>) {
	next if /^#/;
	chomp;
	my @items = split /\t/;
	next unless exists $scaffolds{$items[0]};
	$ScaffoldLen{$items[0]} = $items[1];
}
close I;
print "\n",scalar keys(%ScaffoldLen)," scaffold(s) Load.\nLength:\n";
for (@Scaffolds) {
	my $len = 0;
	for (@$_) {
		$ScaffoldsOffect{$_} = $len;
		$len += $ScaffoldLen{$_};
	}
	#push @t,$len;	# not necessary
	push @LineLen,$len;
	print "[",join(', ',@$_),"]: $len\n";
	$maxLen = $len if $maxLen < $len;
}
print "MaxLen: $maxLen\n";
#ddx \@LineLen;
ddx [\%ScaffoldsOffect,\%ScaffoldLen];
# plothap.pl:65: [[0, 5658748, 12958355], [0, 19124984]]
# plothap.pl:63: { scaffold1458 => 5658748, scaffold188 => 0, scaffold75 => 0 }

my %LD;
open I,'<',$inf or die $!;
while (<I>) {
	next if /^Locus1\t/;
	chomp;
	my @items = split /\t/;
	die scalar @items if @items != 17;
	my ($chr1,$pos1,$id1,undef,undef,undef,$chr2,$pos2,$id2,undef,undef,undef,$dist,$rsq,$dp,$p,$N) = @items;
	my $x = [$chr1,$pos1];
	my $y = [$chr2,$pos2];
#ddx [$x,$y];
	($x,$y) = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } ($x,$y);
#ddx [$x,$y];
	$LD{$$x[0]}{$$x[1]}{$$y[0]}{$$y[1]} = [$dist,$rsq,$dp,$p,$N];
}
#ddx \%LD;

my $Xrange = 500;
my $Yrange = $Xrange/2 + 50;
my $YmaxVal = 2;
my $ArrowLen = 20;	# 16+4
my $axisTick = 4;
my $OutBorder = 24;
my $InBorder = 40;
my $Xtotal = $Xrange + 2*$OutBorder;
my $Yitem = $Yrange + $InBorder;
my $FontSize = int($Xrange/40);
my $FontFamily = 'Arial';

my $perUnit = int($maxLen/10);	# 279.330936 M /10 = 27.933093 M
my $numlevel = int(log($perUnit)/log(10));	# 7
my $numSuflevel = int($numlevel/3);	# 2
my $numSuf=( '', qw( K M G T P E Z Y ) )[$numSuflevel];	# M <---
$numSuflevel = 10 ** (3*$numSuflevel);	# 1,000,000 <---
my $roundTo = 5 * (10 ** ($numlevel-1));	# 5e6
my $unit = $perUnit + (-$perUnit % $roundTo);	# 30 M
my $countMax = int($maxLen/$unit) + (($maxLen%$unit)?1:0);
my $BasepPx = 10*$unit/$Xrange;
#print join(",",$BasepPx,$maxLen/$numSuflevel,$perUnit/$numSuflevel,$numlevel,$numSuf,$numSuflevel,$roundTo/$numSuflevel,$unit/$numSuflevel,$countMax),"\n";
my $ChrCount = @Scaffolds;
my $Ytotal = $Yitem*$ChrCount - $InBorder + 2*$OutBorder;

sub getVal($) {
	my @dat = @{$_[0]};
	my ($sum,$cnt,$max)=(0,0,0);
	for my $k (@dat) {
		$sum += $k;
		++$cnt;
	}
	if ($cnt) {
		return $sum/$cnt;
	} else {
		return -1;
	}
}

my $GRIDsize = 5;
my %PlotLD;
for (@Scaffolds) {
	my %locushere = map { $scaffolds{$_} => 1 } @$_;
	for my $scaff (@$_) {
		my $locus1 = $scaffolds{$scaff};
		for my $pos1 (sort {$a<=>$b} keys %{$LD{$locus1}}) {
			for my $locus2 (keys %{$LD{$locus1}{$pos1}}) {
				next unless exists $locushere{$locus2};
				for my $pos2 (sort {$a<=>$b} keys %{$LD{$locus1}{$pos1}{$locus2}}) {
					my @Offects = map {$ScaffoldsOffect{$scaff2chr{$_}} } ($locus1,$locus2);
#print join(',',$scaff,$locus1,$pos1,$locus2,$pos2,@Offects),"\n";
					my $bp1 = $pos1 + $Offects[0];
					my $bp2 = $pos2 + $Offects[1];
					my ($p1,$p2) = map { int($_/($BasepPx*$GRIDsize))*$GRIDsize } ($bp1,$bp2);
#print join(',',$scaff,$bp1,$bp2,$p1,$p2),"\n";
					my $val = $LD{$locus1}{$pos1}{$locus2}{$pos2}->[1];	# r^2
					next if $val eq "NaN";
					push @{$PlotLD{$locus1}{$p1}{$p2}},$val;
				}
			}
		}
	}
}
#ddx \%PlotLD;

open O,'>',$inf.$outf.'.svg' or die $!;
print O <<HEAD;
<?xml version="1.0"?>
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.2" baseProfile="tiny"
 width="${Xtotal}pt" height="${Ytotal}pt" viewBox="0 0 $Xtotal $Ytotal">
<title>Plot $inf</title>
  <rect x="0" y="0" width="$Xtotal" height="$Ytotal" fill="none" stroke="red" stroke-width="2" />
  <defs>
HEAD

$t=0;
for (@LineLen) {
	my $len = $_ / $BasepPx;
	print O <<CLIP;
    <clipPath id="curveClip$t">
      <path id="curve$t" d="M$len 0 H0 V$len z"/>
    </clipPath>
CLIP
	++$t;
}
print O "  </defs>\n";
$t=0;
my $theY=$OutBorder;
for (@LineLen) {
	my $len = $_ / $BasepPx;
	my $intlen = int(0.999 + $len);
	my $halflen = $intlen/2;
	my $Ylen = $halflen + 50;
	print O <<DEF2;
  <g transform="translate($OutBorder,$theY)" stroke-width="2" stroke="black" font-size="$FontSize" font-family="$FontFamily">
  <rect x="0" y="0" width="$Xrange" height="$Ylen" fill="none" stroke="navy" stroke-width="2" />
    <g transform="matrix(0.5,0.5,-0.5,0.5,$halflen,0)" clip-path="url(#curveClip$t)" stroke-width="0">
      <rect x="0" y="0" width="$len" height="$len" fill="grey" stroke="green" stroke-width="1" />
DEF2
	;
	for my $scaff (@{$Scaffolds[$t]}) {
		my $locus1 = $scaffolds{$scaff};
		for my $pos1 (sort {$a<=>$b} keys %{$PlotLD{$locus1}}) {
			for my $pos2 (sort {$a<=>$b} keys %{$PlotLD{$locus1}{$pos1}}) {
				my $val = getVal($PlotLD{$locus1}{$pos1}{$pos2});
				next if $val == -1;
				my $color = colormap($val);
#print join(',',$locus1,$pos1,$pos2,$val),"\n";
				my $py = $intlen - $pos2 - $GRIDsize;
				print O <<GRID;
      <rect x="$pos1" y="$py" width="$GRIDsize" height="$GRIDsize" fill="$color"/>
GRID
			}
		}
	}

	print O '    </g>
      <g transform="translate(0,',5 + $len/2,')">
      <line x1="0" y1="0" x2="',$len,'" y2="0" stroke-width="3"/>
';
	unless ($t) {
		print O '      <line x1="0" y1="0" x2="',$ScaffoldLen{'scaffold75'}/$BasepPx,'" y2="0" stroke="green" stroke-width="3"/>
';
	}
	print O "      </g>\n  </g>\n";
	++$t;
	$theY += $Ylen + $InBorder;
}
print O "</svg>\n";
close O;

ddx \%Stat;

__END__

./plothap.pl ../rec.npa 17all75.txt

