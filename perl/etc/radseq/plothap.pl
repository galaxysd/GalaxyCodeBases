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

my ($maxLen,@LineLen,%ScaffoldLen)=(0);
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
	$len += $ScaffoldLen{$_} for @$_;
	push @LineLen,$len;
	print "[",join(', ',@$_),"]: $len\n";
	$maxLen = $len if $maxLen < $len;
}
print "MaxLen: $maxLen\n";
#ddx \@LineLen;

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
my $Yrange = int($Xrange*sqrt(2)) + 120;
my $YmaxVal = 5;
my $ArrowLen = 20;	# 16+4
my $axisTick = 4;
my $OutBorder = 24;
my $InBorder = 40;
my $Xtotal = $Xrange + $ArrowLen + 2*$OutBorder;
my $Yitem = $Yrange + $ArrowLen + $InBorder;
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

sub getVal($) {
	my %dat = %{$_[0]};
	my ($sum,$cnt,$max,$major)=(0,0,0,0);
	for my $k (keys %dat) {
		$sum += $k * $dat{$k};
		$cnt += $dat{$k};
	}
	if ($cnt) {
		return $sum/$cnt;
	} else {
		return -1;
	}
}

my %PlotLD;
for (@Scaffolds) {
	for my $scaff (@$_) {
		my $locus1 = $scaffolds{$scaff};
		for my $pos1 (sort {$a<=>$b} keys %{$LD{$locus1}}) {
			for my $locus2 (keys %{$LD{$locus1}{$pos1}}) {
				for my $pos2 (sort {$a<=>$b} keys %{$LD{$locus1}{$pos1}{$locus2}}) {
					#print join(',',$scaff,$locus1,$pos1,$locus2,$pos2),"\n";
					my ($p1,$p2) = map { int(10*$_/$BasepPx)/10 } ($pos1,$pos2);
					my $val = $LD{$locus1}{$pos1}{$locus2}{$pos2}->[1];
					push @{$PlotLD{$locus1}{$p1}{$locus2}{$p2}},$val;
				}
			}
		}
	}
}
ddx \%PlotLD;

__END__

./plothap.pl ../rec.npa 17all75.txt

