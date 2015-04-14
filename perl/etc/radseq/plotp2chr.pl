#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Data::Dump qw(ddx);
use Galaxy;

die "Usage: $0 <prefix> [out midfix]\n" if @ARGV<1;
my $inf=shift;
my $outf=shift;

if ($outf) {
	$outf = ".$outf";
} else {
	$outf = '';
}

my $scaffnfo = 'chr.nfo';
$scaffnfo = 'chr.inf';
my $markerdat = '/share/users/huxs/work/tiger/paper/rec.npa';
$markerdat = "/share/users/huxs/work/tiger/cat/rec.npa";
$markerdat = $inf . '.npa';
print "From [$inf] to [$inf$outf.(dat|svg)]\n";

my %MVScaffolds = map {$_ => 1} qw(scaffold75 scaffold1458 scaffold188);

my %Stat;

sub getVal($) {	# deprecated
	my %dat = %{$_[0]};
	my ($sum,$cnt,$max,$major)=(0,0,0,0);
	my $mjcnt=0;
	for my $k (keys %dat) {
		$sum += $k * $dat{$k};
		$cnt += $dat{$k};
		$max = $k if $max < $k;
		if ($k and $mjcnt<$dat{$k}) {
			$major = $k;
			$mjcnt = $dat{$k};
		}
	}
	if ($cnt) {
		return ($major,$max);
		return $sum/$cnt;
	} else {
		return -1;
	}
}
sub getCircles($) {
	my %dat = %{$_[0]};
	my @ret;
	my ($max,$maxR) = (0,0);
	for my $k (keys %dat) {
		next unless $k;
		my $r = int( 100 * sqrt($dat{$k}) ) / 100;
		$r = 1.35 if $r > 1.35;
		push @ret,[$k,$r];
		if ($max < $k) {
			$max = $k;
			$maxR = $dat{$k};
		}
	}
	unshift @ret,[$max,$maxR];
	return @ret;
}

my (%ScaffoldLen,@ChrOrder0,%ChrID,%ChrColor);
open I,'<',$scaffnfo or die $!;
while (<I>) {
	next if /^#/;
	chomp;
	my @items = split /\t/;
	$ScaffoldLen{$items[0]} = $items[1];
	$ChrID{$items[0]} = $items[2];
	$ChrColor{$items[0]} = $items[3];
	push @ChrOrder0,$items[0];
}
close I;
print "\n",scalar keys(%ScaffoldLen)," scaffold(s) Load.\n";

my %MarkerDat;
open I,'<',$markerdat or die $!;
while (<I>) {
	next if /^#/;
	chomp;
	my @items = split /\t/;
	die unless exists $ScaffoldLen{$items[1]};
    my @tmp = split /[\/|]/,$items[7].'/'.$items[8];    # Aff,UnAff
	push @{$MarkerDat{$items[1]}},[ $items[2],$items[9],int(0.5-1000*log($items[9])/log(10))/1000,$tmp[1]+$tmp[2] ];	# pos,p,lg(p) round to 0.001,SumMid(0 is red)
}
close I;
#ddx \%MarkerDat;
my (%OrderbyChr,%OrderedOnce,@DirectOrder,@notOrdered);
open I,'<',$inf.'.order' or warn $! and goto NULLORDERFILE;
while (<I>) {
	next if /^#/;
	chomp;
	my ($scaff,$chr) = split /\t/ or next;	# skip blank lines
	next if exists $OrderedOnce{$scaff};
	push @{$OrderbyChr{$chr}},$scaff;
	push @DirectOrder,$scaff;
	++$OrderedOnce{$scaff};
}
close I;
NULLORDERFILE:
my $TotalLen=0;
#my @order = map {$_='scaffold'.$_} sort {$a<=>$b} map {s/^scaffold//;$_;} keys %MarkerDat;
my @order = @ChrOrder0;
for my $scaff (@order) {
	next unless exists $MarkerDat{$scaff};
	if (exists $OrderedOnce{$scaff}) {
		$Stat{Len_Ordered} += $ScaffoldLen{$scaff} or die;
		++$Stat{Scaffold_Ordered};
	} else {
		push @notOrdered,$scaff;
		$Stat{Len_notOrdered} += $ScaffoldLen{$scaff} or die;
		++$Stat{Scaffold_notOrdered};
	}
	$TotalLen += $ScaffoldLen{$scaff};
}
if ($Stat{Scaffold_notOrdered}) {
	$Stat{AvgLen_notOrdered} = $Stat{Len_notOrdered} / $Stat{Scaffold_notOrdered};
}
if ($Stat{Scaffold_Ordered}) {
	$Stat{AvgLen_Ordered} = $Stat{Len_Ordered} / $Stat{Scaffold_Ordered};
}

# ------ BEGIN PLOT --------
# 1in = 2.54cm = 25.4mm = 72pt = 12pc, 1pc=2.1167mm, 1pt=0.35278mm
my @color = qw(Red Purple Brown Navy Green Maroon Blue Teal);
my $Xrange = 500;
my $Yrange = 320;
my $YmaxVal = 5;
my $ArrowLen = 20;	# 16+4
my $axisTick = 4;
my $OutBorder = 24;
my $InBorder = 40;
my $Yextra = 200;
my $Xtotal = $Xrange + $ArrowLen + 2*$OutBorder;
my $Yitem = $Yrange + $ArrowLen + $InBorder + $Yextra;
my $FontSize = int($Xrange/40);
my $FontFamily = 'Arial';

my $perUnit = int($TotalLen/10);	# 279.330936 M /10 = 27.933093 M
my $numlevel = int(log($perUnit)/log(10));	# 7
my $numSuflevel = int($numlevel/3);	# 2
my $numSuf=( '', qw( K M G T P E Z Y ) )[$numSuflevel];	# M <---
$numSuflevel = 10 ** (3*$numSuflevel);	# 1,000,000 <---
my $roundTo = 5 * (10 ** ($numlevel-1));	# 5e6
my $unit = $perUnit + (-$perUnit % $roundTo);	# 30 M
my $countMax = int($TotalLen/$unit) + (($TotalLen%$unit)?1:0);
#print join(",",$TotalLen/$numSuflevel,$perUnit/$numSuflevel,$numlevel,$numSuf,$numSuflevel,$roundTo/$numSuflevel,$unit/$numSuflevel,$countMax),"\n";
my $BasepPx = 10*$unit/$Xrange;

my @Yticks;
for my $i (1 .. $YmaxVal) {	# 0 will be shared with X-axie
	my $pY = $Yrange - $i * ($Yrange/$YmaxVal);
	push @Yticks,$pY;
}

my (%PlotDat,%PlotScaffRange,%PlotCandiDat);
my $start = 0;
for my $scaff (@notOrdered,@DirectOrder) {	# Well, Ordered last to be far away from X-axie.
	unless (exists $MarkerDat{$scaff}) {
		warn "[!marker]$scaff\n";
		++$Stat{'Marker_Not_found'};
		next;
	}
	my $maxlgp = 0;
	for my $mditem (@{$MarkerDat{$scaff}}) {
		my ($pos,$p,$lgp,$isFalse) = @$mditem;
		my $posOchr = $start + $pos;
		if ($posOchr < 0) {
			$posOchr = 0;
			++$Stat{'Marker_Pos_Minus'};
		} elsif ($pos > $ScaffoldLen{$scaff}) {
			$posOchr = $start + $ScaffoldLen{$scaff};
			++$Stat{'Marker_Pos_Overflow'};
		}
        my $value = int(0.5+10*$posOchr/$BasepPx)/10;	# 10 => 720 dpi for pt unit system, enough.
		++$PlotDat{$scaff}{$value}{$lgp};
        ++$PlotCandiDat{$scaff}{$value}{$lgp} unless $isFalse;
		$maxlgp = $lgp if $maxlgp < $lgp;
	}
	$PlotScaffRange{$scaff} = [ ( map {int(0.5+10*$_/$BasepPx)/10} ($start+1,$start+$ScaffoldLen{$scaff}) ),$maxlgp ];
	$start += $ScaffoldLen{$scaff};
}

my $ChrCount = keys %PlotDat;
$ChrCount = 1;
my $Ytotal = $Yitem*$ChrCount - $InBorder + $ArrowLen + 2*$OutBorder;

open O,'>',$inf.$outf.'.svg' or die $!;
print O <<HEAD;
<?xml version="1.0"?>
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.2" baseProfile="tiny"
 width="${Xtotal}pt" height="${Ytotal}pt" viewBox="0 0 $Xtotal $Ytotal">
<title>Plot $inf</title>
<!--
  <rect x="0" y="0" width="$Xtotal" height="$Ytotal" fill="none" stroke="red" stroke-width="2" />
-->
HEAD

print O <<DEF1;
  <defs>
    <g id="axis" stroke="black" font-size="$FontSize" font-family="$FontFamily" stroke-width="0" text-anchor="middle">
      <polyline fill="none" points="-2,-4 0,-20 2,-4 0,-20 
DEF1
for (@Yticks) {
	print O "        0,$_ -$axisTick,$_ 0,$_ \n"
}
print O "        0,$Yrange\n";

for my $i (1 .. 10) {
	my $x = $i*$Xrange/10;
	print O '        ',$x,',',$Yrange,' ',$x,',',$Yrange+$axisTick,' ',$x,',',$Yrange,"\n";
}
print O '        ',$Xrange+$ArrowLen,',',$Yrange,' ',
	$Xrange+$axisTick,',',$Yrange-2,' ',$Xrange+$axisTick,',',$Yrange+2,' ',
	$Xrange+$ArrowLen,',',$Yrange,'" stroke-width="2"/>',"\n";
for my $i (1 .. $YmaxVal) {
	my $y = $Yticks[$i-1];
	print O <<TXTAX;
      <text x="0" y="$y" dx="-20" dy="5" text-anchor="middle">$i</text>
      <line x1="0" y1="$y" x2="$Xrange" y2="$y" stroke-width="0.5" stroke-dasharray="5,9"/>
TXTAX
}
for my $i (0 .. 10) {
	my $x = $i*$Xrange/10;
	my $l = $unit*$i/$numSuflevel;
	$l .= " $numSuf" if $l;
	print O <<TXT1;
      <text x="$x" y="$Yrange" dy="20" text-anchor="middle">$l</text>
TXT1
}
print O <<DEF2;
    </g>
  </defs>
  <g transform="translate($OutBorder,$OutBorder)" stroke-width="2" stroke="black" font-size="$FontSize" font-family="$FontFamily">
DEF2

my %maxCircles;
my $thisChrNo=0;
for my $chr ('-lg(p)') {
	my $topY = $ArrowLen + $Yitem*$thisChrNo;
	print O <<TXT2;
    <g transform="translate(0,$topY)" stroke-width="1">
      <use x="0" y="0" xlink:href="#axis" stroke-width="2"/>
      <text x="8" y="-10" stroke-width="0">$chr</text>
TXT2
	my $scaffcnt=0;
	for my $scaff (@DirectOrder,@notOrdered) {
		next unless exists $MarkerDat{$scaff};
		my @Poses = sort {$a<=>$b} keys %{$PlotDat{$scaff}};
		my $thiscolor = $color[$scaffcnt%scalar(@color)];
		my ($pa,$pb,$maxlgp) = @{$PlotScaffRange{$scaff}};
# 		if (exists $MVScaffolds{$scaff}) {
# 			$thiscolor = 'red';
# 		} else {
# 			$thiscolor = 'navy';
# 		}
		$thiscolor = '#'.$ChrColor{$scaff};
		print O '      <g stroke="',$thiscolor,'" fill="',$thiscolor,'" focusable = "true">',"\n        <title>$scaff, max=$maxlgp</title>\n";
		my $topestY = $Ytotal;
		for my $pos (@Poses) {
			my @Circles = getCircles($PlotDat{$scaff}{$pos});
            my @CandiCircles = getCircles($PlotCandiDat{$scaff}{$pos}) if exists $PlotCandiDat{$scaff}{$pos};
			my $t = shift @Circles;
            shift @CandiCircles;
			my $p = int($pos);
			unless (exists $maxCircles{$p}) {
				$maxCircles{$p} = $t;
			} else {
				if ($maxCircles{$p}->[0] < $t->[0]) {
					$maxCircles{$p} = $t;
				} elsif ($maxCircles{$p}->[0] == $t->[0] and $maxCircles{$p}->[1] < $t->[1]) {
					$t->[1] += $maxCircles{$p}->[1];
					$maxCircles{$p} = $t;
				}
			}
			for (@Circles) {
				my ($y,$r) = @$_;
				my $Py = int(10*$Yrange*(1-$y/$YmaxVal))/10;
				$topestY = $Py if $topestY > $Py;
#print "$y,$Py,$Yrange,$YmaxVal\n";
				print O "        <circle cx=\"$pos\" cy=\"$Py\" r=\"$r\" />\n";
			}
			for (@CandiCircles) {
				my ($y,$r) = @$_;
				my $Py = int(10*$Yrange*(1-$y/$YmaxVal))/10;
				$topestY = $Py if $topestY > $Py;
                #print "$y,$Py,$Yrange,$YmaxVal\n";
				print O "        <circle cx=\"$pos\" cy=\"$Py\" r=\"$r\" stroke=\"red\" fill=\"red\" />\n";
			}
=pod
			my ($pYmajor,$pYmax) = map {int(10*$Yrange*(1-$_/$YmaxVal))/10;} ($major,$max);
			print O <<TXTL;
        <circle cx="$pos" cy="$pYmajor" r="1" />
        <circle cx="$pos" cy="$pYmax" r="0.5" />
TXTL
=cut
		}
		#my ($pa,$pb) = @{$PlotScaffRange{$scaff}};
		print O <<TXTLB;
        <line x1="$pa" y1="$Yrange" x2="$pb" y2="$Yrange" stroke-width="2"/>
      </g>
TXTLB
		print O "      <text x=\"$Poses[0]\" y=\"",$topestY-12,"\" dy=\"5\" text-anchor=\"middle\" fill=\"navy\" stroke-width=\"0\">$scaff</text>\n" if exists $MVScaffolds{$scaff};
		++$scaffcnt;
	}
#	print O '      <polyline fill="none" stroke="gold" stroke-width="0.5" points="';
#	for my $x (sort {$a<=>$b} keys %maxCircles) {
#		my ($y,$s) = @{$maxCircles{$x}};
#		next unless $y;
#		$y = int(10*$Yrange*(1-$y/$YmaxVal))/10;
#		print O "$x,$y ";
#	}
#	print O "\" />\n    </g>\n";
	my ($xx,$yy)=(0,$Yrange + 30);
	for my $id (@ChrOrder0) {
		my $name = $ChrID{$id};
		my $color = '#'.$ChrColor{$id};
		my $tx = $xx + 17;
		print O <<CHRBOX;
	<rect x="$xx" y="$yy" rx="5" ry="5" width="15" height="15" style="fill:$color;stroke:black;stroke-width:1"/>
	<text x="$tx" y="$yy" dy="12" stroke-width="0">Chr$name</text>
CHRBOX
	$xx += 56;
	if ($xx > $Xtotal - 70) {
		$yy += 24;
		$xx = 0;
	}
	}

	print O "    </g>\n";
#	++$thisChrNo;
}

print O "  </g>\n</svg>\n";
=pod
print O "
  </g>
  <rect x=\"$OutBorder\" y=\"$OutBorder\" width=\"",$Xrange+$ArrowLen,"\" height=\"",$Ytotal-2*$OutBorder,"\" fill=\"none\" stroke=\"blue\" stroke-width=\"1\" />
  <line x1=\"",$Xrange+$OutBorder,"\" y1=\"$OutBorder\" x2=\"",$Xrange+$OutBorder,"\" y2=\"",$Ytotal-$OutBorder,"\" stroke=\"blue\" stroke-width=\"1\"/>
</svg>
";
=cut
close O;

ddx \%Stat;
print commify($TotalLen),"\n";
__END__
perl -lane '$F[-5]=~s/\,$//;print join("\t",@F[0,1,-5],'f00')' /bak/seqdata/genomes/Felis_catus-6.2/chr.nfo > chr.inf

perl -F',' -lane "map {s/'//g;s/ /_/} @F;print join(\"\t\",@F);" ./tig2cat.csv > ./tig2cat.tsv
sort -t $'\t' -k 1,1 -k 2,2n -k 3,3n ./tig2cat.tsv > ./tig2cat.tsv.s
awk '{print "scaffold"$4"\tchr"$1}' tig2cat.tsv.s|uniq > tig2cat.order

./plotscaf.pl rec16q20.npa15 tig2cats16q20s15
