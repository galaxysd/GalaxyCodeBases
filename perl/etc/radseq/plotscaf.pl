#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Data::Dump qw(ddx);
use Galaxy;

die "Usage: $0 <order file> <marker.rec> <outfile> <TRUE for 97&1457 only,also dot size>\n" if @ARGV<3;
my $orderf=shift;
my $markerdat=shift;
my $outf=shift;
my $type=shift;

my $scaffnfo = '/bak/seqdata/genomes/TigerRefGenome/chr.nfo';
#$scaffnfo = '/bak/seqdata/genomes/Felis_catus-6.2/chr.nfo';

print "From [O:$orderf][$markerdat] to [$outf.(dat|svg)]\n";

my %MVScaffolds = map {$_ => 1} qw(scaffold97 scaffold1457 scaffold91);
if ($type) {
	%MVScaffolds = ();
}
my %ResultZone = ('scaffold97' => -1, 'scaffold1457' => 1);	# scaffold97 on '-' strand.
my $CHR_UN = 'UN';

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

my %CirclesR;
my $MaxCirclesR = 1.5;
sub getCircles($) {
	my %dat = %{$_[0]};
	my @ret;
	my ($max,$maxR) = (0,0);
	for my $k (keys %dat) {
		next unless $k;
		my $r = int( 100 * sqrt($dat{$k}) ) / 100;
		++$CirclesR{$r};
		$r = $MaxCirclesR if $r > $MaxCirclesR;
		push @ret,[$k,$r];
		if ($max < $k) {
			$max = $k;
			$maxR = $dat{$k};
		}
	}
	unshift @ret,[$max,$maxR];
	return @ret;
}

my (%ScaffoldLen,@ChrOrder0);
open I,'<',$scaffnfo or die $!;
while (<I>) {
	next if /^#/;
	chomp;
	my @items = split /\t/;
	$ScaffoldLen{$items[0]} = $items[1];
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
	if ($type) {
		next unless $ResultZone{$items[1]};
	}
	die unless exists $ScaffoldLen{$items[1]};
    my @tmp = split /[\/|]/,$items[7].'/'.$items[8];    # Aff,UnAff
	push @{$MarkerDat{$items[1]}},[ $items[2],$items[9],int(0.5-1000*log($items[9])/log(10))/1000,$tmp[1]+$tmp[2] ];	# pos,p,lg(p) round to 0.001,SumMid(0 is red)
}
close I;
# ddx \%MarkerDat;
#   scaffold1457 => [
#                     [8766, 0.3698, 0.432, 7],
#                     [33512, 0.02941, 1.532, 4],
#                     [34634, 0.4706, 0.327, 7],
#                     [36213, 0.03251, 1.488, 5],
#                     [6453045, 1, 0, 9],
#                   ],
my (%OrderbyChr,@ChrOrder,%OrderedOnce,@DirectOrder,@notOrdered);
my (@order,%Scaff2Chr);
open I,'<',$orderf or warn $! and goto NULLORDERFILE;
while (<I>) {
	next if /^#/;
	chomp;
	my ($scaff,$chr) = split /\t/ or next;	# skip blank lines
	next if exists $OrderedOnce{$scaff};
	next if $type and !$ResultZone{$scaff};	# only those in %ResultZone
	push @ChrOrder,$chr unless exists $OrderbyChr{$chr};
	push @{$OrderbyChr{$chr}},$scaff;
	push @DirectOrder,$scaff;
	++$OrderedOnce{$scaff};
	$Scaff2Chr{$scaff} = $chr;
}
close I;

@order = @ChrOrder0;
### my %UsedChr = map {$_ => 0} @ChrOrder;
goto TOCONTINUE;

NULLORDERFILE:
@order = map {$_='scaffold'.$_} sort {$a<=>$b} map {s/^scaffold//;$_;} keys %MarkerDat;

TOCONTINUE:
my $TotalLen=0;
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
	### $UsedChr{ $Scaff2Chr{$scaff} } = 1 if $Scaff2Chr{$scaff};	# Well, assume the %ResultZone is mapped.
}
if ($Stat{Scaffold_notOrdered}) {
	$Stat{AvgLen_notOrdered} = $Stat{Len_notOrdered} / $Stat{Scaffold_notOrdered};
}
if ($Stat{Scaffold_Ordered}) {
	$Stat{AvgLen_Ordered} = $Stat{Len_Ordered} / $Stat{Scaffold_Ordered};
}
=pod
if ($type) {
	my @chrs;
	for (@ChrOrder) {
		push @chrs,$_ if $UsedChr{$_};
	}
	@ChrOrder = @chrs;
	ddx \@notOrdered;
}
=cut
push @ChrOrder,$CHR_UN unless $type;

# ------ BEGIN PLOT --------
# 1in = 2.54cm = 25.4mm = 72pt = 12pc, 1pc=2.1167mm, 1pt=0.35278mm
my @color = qw(Black Red Green Navy Blue Purple Orange Gray Maroon Teal Brown);
@color = qw(Black Brown #0F8B43 #3954A5 #199BCD #B2499B #EE7821 #686868);	# #ED1924
my $Xrange = 1000;
my $Yrange = 320;
$Xrange = 360 if $type;
my $YmaxVal = 6;
my $ArrowLen = 20;	# 16+4
my $axisTick = 4;
my $OutBorder = 24;
my $InBorder = 40;
my $Xtotal = $Xrange + $ArrowLen + 2*$OutBorder;
my $Yitem = $Yrange + $ArrowLen + $InBorder;
my $FontSize = int($Yrange/16);
my $SmallFontSize = $FontSize-2;
my $FontFamily = 'Arial';

my $perUnit = int($TotalLen/10);	# 279.330936 M /10 = 27.933093 M
my $numlevel = int(log($perUnit)/log(10));	# 7
my $numSuflevel = int($numlevel/3);	# 2
my $numSuf=( '', qw( K M G T P E Z Y ) )[$numSuflevel];	# M <---
$numSuflevel = 10 ** (3*$numSuflevel);	# 1,000,000 <---
my $roundTo = 5 * (10 ** ($numlevel-1));	# 5e6
my $unit = $perUnit + (-$perUnit % $roundTo);	# 30 M
my $countMax = int($TotalLen/$unit) + (($TotalLen%$unit)?1:0);
print "[!]PlotSize: ",join(",",$TotalLen/$numSuflevel,$perUnit/$numSuflevel,$numlevel,$numSuf,$numSuflevel,$roundTo/$numSuflevel,$unit/$numSuflevel,$countMax),"\n";
my $BasepPx = 10*$unit/$Xrange;
my $Ymin = 0;

my @Yticks;
for my $i ($Ymin .. $YmaxVal) {	# 0 will be shared with X-axie, set $Ymin=1 if needed this.
	my $pY = $Yrange - $i * ($Yrange/$YmaxVal);
	push @Yticks,$pY;
}

my (%PlotDat,%PlotScaffRange,%PlotCandiDat);
my $start = 0;
for my $scaff (@DirectOrder,@notOrdered) {	# Well, I prefer Ordered last to be far away from Y-axie, but they does not.
	unless (exists $MarkerDat{$scaff}) {
		warn "[!marker]$scaff\n" unless $type;
		++$Stat{'Marker_Not_found'};
		next;
	}
	my $maxlgp = 0;
	if ($type && $ResultZone{$scaff} < 0) {
		my @tmpDat = @{$MarkerDat{$scaff}};
		my %poses;
		for my $mditem (@tmpDat) {
			my ($pos,$p,$lgp,$isFalse) = @$mditem;
			++$poses{$pos};
		}
		my ($maxpos,$minpos) = (sort {$b<=>$a} keys %poses)[0,-1];
		my $ToDeltra = $maxpos+$minpos;
		$MarkerDat{$scaff} = [];
		for my $mditem (@tmpDat) {
			my ($pos,$p,$lgp,$isFalse) = @$mditem;
			push @{$MarkerDat{$scaff}},[$ToDeltra-$pos,$p,$lgp,$isFalse];
		}
	}
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
        #++$PlotCandiDat{$scaff}{$value}{$lgp} unless $isFalse;
		$maxlgp = $lgp if $maxlgp < $lgp;
	}
	$PlotScaffRange{$scaff} = [ ( map {int(0.5+10*$_/$BasepPx)/10} ($start+1,$start+$ScaffoldLen{$scaff}) ),$maxlgp ];
	$start += $ScaffoldLen{$scaff};
}

my $ChrCount = keys %PlotDat;
$ChrCount = 1;
my $Ytotal = $Yitem*$ChrCount - $InBorder + $ArrowLen + 2*$OutBorder;

open O,'>',$outf.'.svg' or die $!;
print O <<HEAD;
<?xml version="1.0"?>
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.2" baseProfile="tiny"
 width="${Xtotal}pt" height="${Ytotal}pt" viewBox="0 0 $Xtotal $Ytotal">
<title>Plot $markerdat [O:$orderf]</title>
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

my $XaxisColor = 'white';
my $XaxisCount = 10;
if ($type) {
	$XaxisColor = 'black';
	$XaxisCount = 5;
} else {
	$axisTick=1;	# remove X ticks.
}
for my $i (1 .. $XaxisCount) {
	my $x = $i*$Xrange/$XaxisCount;
	print O '        ',$x,',',$Yrange,' ',$x,',',$Yrange+$axisTick,' ',$x,',',$Yrange,"\n";
}
print O '        ',$Xrange+$ArrowLen,',',$Yrange,' ',
	$Xrange+$axisTick,',',$Yrange-2,' ',$Xrange+$axisTick,',',$Yrange+2,' ',
	$Xrange+$ArrowLen,',',$Yrange,'" stroke-width="2"/>',"\n";
for my $i ($Ymin .. $YmaxVal) {
	my $y = $Yticks[$i-$Ymin];
	print O <<TXTAX;
      <text x="0" y="$y" dx="-20" dy="5" text-anchor="middle">$i</text>
      <line x1="0" y1="$y" x2="$Xrange" y2="$y" stroke-width="0.5" stroke-dasharray="5,9" stroke="$XaxisColor"/>
TXTAX
}

for my $i (0 .. $XaxisCount) {
	my $x = $i*$Xrange/$XaxisCount;
	my $l = (10/$XaxisCount) * $unit*$i/$numSuflevel;
	$l .= " $numSuf" if $l;
	print O <<TXT1;
      <text x="$x" y="$Yrange" dy="20" text-anchor="middle" font-size="$SmallFontSize" fill="$XaxisColor">$l</text>
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
	my $chrcnt=0;
	for my $chr ( @ChrOrder ) {
		my ($pChrA,$pChrB)=($Xrange,0);
		my @scaffs;	# @DirectOrder,@notOrdered
		if ($chr ne $CHR_UN) {	# 'UN'
			@scaffs = @{$OrderbyChr{$chr}};
		} else {
			@scaffs = @notOrdered;
		}
		my $thiscolor = $color[$chrcnt%scalar(@color)];
		print O '      <g stroke="',$thiscolor,'" fill="',$thiscolor,'" focusable = "true">',"\n        <title>$chr</title>\n";
		for my $scaff ( @scaffs ) {
			next unless exists $MarkerDat{$scaff};
			my @Poses = sort {$a<=>$b} keys %{$PlotDat{$scaff}};
			#my $thiscolor = $color[$scaffcnt%scalar(@color)];
			my ($pa,$pb,$maxlgp) = @{$PlotScaffRange{$scaff}};
			$pChrA = $pa if $pChrA > $pa;
			$pChrB = $pb if $pChrB < $pb;

			if (exists $MVScaffolds{$scaff}) {
				$thiscolor = 'red';
			} else {
				$thiscolor = $color[$chrcnt%scalar(@color)];
			}
#			$thiscolor = 'black';

			print O '      <g stroke="',$thiscolor,'" fill="',$thiscolor,'" focusable = "true">',"\n        <title>$scaff, max=$maxlgp</title>\n";
			my $topestY = $Ytotal;
			for my $pos (@Poses) {
				my @Circles = getCircles($PlotDat{$scaff}{$pos});
	            #my @CandiCircles = getCircles($PlotCandiDat{$scaff}{$pos}) if exists $PlotCandiDat{$scaff}{$pos};
				my $t = shift @Circles;
	            #shift @CandiCircles;
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
					if ($type) {
						$r *= $type;
					} else {
						$r *= 2;	# make it bigger
					}
					my $Py = int(10*$Yrange*(1-$y/$YmaxVal))/10;
					$topestY = $Py if $topestY > $Py;
	#print "$y,$Py,$Yrange,$YmaxVal\n";
					#print O "        <circle cx=\"$pos\" cy=\"$Py\" r=\"$r\" />\n";
					print O "        <rect x=\"$pos\" y=\"$Py\" width=\"$r\" height=\"$r\" />\n";
				}
=pod
				for (@CandiCircles) {
					my ($y,$r) = @$_;
					my $Py = int(10*$Yrange*(1-$y/$YmaxVal))/10;
					$topestY = $Py if $topestY > $Py;
	                #print "$y,$Py,$Yrange,$YmaxVal\n";
					print O "        <circle cx=\"$pos\" cy=\"$Py\" r=\"$r\" stroke=\"red\" fill=\"red\" />\n";
				}
=cut
			}
			#my ($pa,$pb) = @{$PlotScaffRange{$scaff}};
			print O "      </g>\n";
			if (exists $MVScaffolds{$scaff} or $type) {
				print O "      <text x=\"$Poses[0]\" y=\"",$topestY-$FontSize,"\" dy=\"5\" text-anchor=\"middle\" fill=\"navy\" stroke-width=\"0\">$scaff</text>\n";
			}
			#++$scaffcnt;
			++$chrcnt if $type;
		}
		$chr =~ s/^chr//i;
		print O "      <text x=\"",($pChrA+$pChrB)/2,"\" y=\"",$Yrange+$SmallFontSize+(1+$SmallFontSize)*($chrcnt%2),
			"\" dy=\"5\" text-anchor=\"middle\" stroke-width=\"0\" font-size=\"$SmallFontSize\">$chr</text>\n";	# one long line will break SyntaxHighlight ...
		print O <<TXTLB;
        <line x1="$pChrA" y1="$Yrange" x2="$pChrB" y2="$Yrange" stroke-width="2"/>
      </g>
TXTLB
		++$chrcnt unless $type;
	}
#	print O '      <polyline fill="none" stroke="gold" stroke-width="0.5" points="';
#	for my $x (sort {$a<=>$b} keys %maxCircles) {
#		my ($y,$s) = @{$maxCircles{$x}};
#		next unless $y;
#		$y = int(10*$Yrange*(1-$y/$YmaxVal))/10;
#		print O "$x,$y ";
#	}
#	print O "\" />\n    </g>\n";
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
ddx \%CirclesR;
__END__

perl -F',' -lane "map {s/'//g;s/ /_/} @F;print join(\"\t\",@F);" ./tig2cat.csv | awk '{if(NR>1)print}' > ./tig2cat.tsv
sort -t $'\t' -k 1,1 -k 2,2n -k 3,3n ./tig2cat.tsv > ./tig2cat.tsv.s
awk '{print "scaffold"$4"\tchr"$1}' tig2cat.tsv.s|uniq > tig2cat.order

./plotscaf.pl rec16q20.npa15 tig2cats16q20s15

======

zcat tig2cat.csv.gz | perl -F',' -lane "map {s/'//g;s/ /_/} @F;print join(\"\t\",@F);" | awk '{if(NR>1)print}' > tig2cat.tsv
sort -t $'\t' -k 1,1 -k 2,2n -k 3,3n ./tig2cat.tsv > ./tig2cat.tsv.s
awk '{print "scaffold"$4"\tchr"$1}' tig2cat.tsv.s|uniq > tig2cat.order

perl plotscaf.pl tig2cat.order sw000-17R.REC.out.rec.npa sw000
./plotscaf.pl tig2cat.order sw210-17R.REC.out.rec.npa sw210ts 1
