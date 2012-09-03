#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Data::Dump qw(ddx);
use Statistics::LineFit;
use Galaxy::IO;
use Galaxy::SeqTools;

die "Usage: $0 <in> <out.prefix>\n" if @ARGV<2;
my $inf=shift;
my $outf=shift;

my $lineFit = Statistics::LineFit->new();
my (%dat,$t,%out);
open I,'<',$inf or die $!;
while (<I>) {
	chomp;
	my @items = split /\s+/;
	# 0NMid 1chr 2+- 3start 4end 5gene 6incmpl5 7incmpl3 8scaffold 9start 10end 11+- 12gene(same)
	if ( $items[2] eq $items[11] ) {
		$t = 0;
	} else {
		$t = 1;
	}
	push @{$dat{$items[8]}{$items[1]}}, [$t,$items[5],$items[3],$items[4],$items[9],$items[10]];
}
close I;
#ddx \%dat;

open OL,'>',$outf.'.plog' or  die $!;
my %OutDat;
for my $scafd (sort keys %dat) {
	my (%scoreSame,%scoreAbs,%thePos);
	for my $chr (keys %{$dat{$scafd}}) {
		my $alignmtsA = $dat{$scafd}{$chr};
		my (%cntStrand);
		@$alignmtsA = sort {$a->[2] <=> $b->[2]} @$alignmtsA;
		for my $item (@$alignmtsA) {
print OL join("\t",$scafd,$chr,@$item);
			++$cntStrand{$$item[0]};
			my $lenChr = $$item[3] - $$item[2];
			my $lenScafd = $$item[5] - $$item[4];
			my ($a,$b) = sort {$a<=>$b} ($lenChr,$lenScafd);
			if ($a) {
				$t = $a / $b;
			} else {
#print "\t***\n";
				#next;
				$t = 0.1;
			}
			push @{$thePos{$chr}},[$$item[0],$$item[2]+$lenChr/2,$$item[4]+$lenScafd/2];
print OL join("\t",'',$lenScafd,$lenChr,$t),"\n";
			if ($$item[0]) {
				$scoreSame{$chr} -= $t;
			} else {
				$scoreSame{$chr} += $t;
			}
			$scoreAbs{$chr} += $t;
		}
		print OL join(",",$scoreAbs{$chr},$scoreSame{$chr},%cntStrand),"\t-----\n";
	}
	my ($theChr) = sort { $scoreAbs{$b} <=> $scoreAbs{$a} } keys %scoreAbs;
	my $diffStrand;
	if ($scoreSame{$theChr} > 0) {
		$diffStrand = 0;
	} else {
		$diffStrand = 1;
	}
	my (@x,@y);
	for my $item (@{$thePos{$theChr}}) {
		next unless $item->[0] == $diffStrand;
		push @y,$item->[1];	# ChrPos
		push @x,$item->[2];	# ScaffPos
	}
	my $scaffOnChr;
	my ($intercept, $slope, $r2,$sigma,@residuals) = (0,0,-1,-1);
	if (@x > 1) {
		$lineFit->setData(\@x, \@y) or die "Invalid regression data\n";
		($intercept, $slope) = $lineFit->coefficients();
		defined $intercept or die "Can't fit line if x values are all equal";
		$r2 = $lineFit->rSquared();
		$scaffOnChr = $intercept + $slope;
		$sigma = $lineFit->sigma();
		@residuals = $lineFit->residuals();
	} elsif (@x == 1) {
		$scaffOnChr = 1 + $y[0] - $x[0];
	} else {
		die;
	}
	$scaffOnChr = int($scaffOnChr);
	print OL '-' x 10, join(',',
	$theChr,$scaffOnChr,'r2',$r2, 's',$sigma,'k',$slope,'b',$intercept, 'dy',(map {int($_)} @residuals),
	"\n",'-' x 10,
	'x',(map {int($_)} @x),'y',(map {int($_)} @y) ),"\n";
	push @{$OutDat{$theChr}},[$scaffOnChr,$diffStrand,$scafd,$r2,$slope,$sigma];
	#print O join("\t",$theChr,$scaffOnChr,$diffStrand,$scafd,$r2,$slope,$sigma),"\n";
}
close OL;

open O,'>',$outf.'.plst' or  die $!;
print O join("\t",'#ChrID','Pos','isDiffStrand','scaffold','r^2','slope(k)','sigma'),"\n";
for my $chr (sort keys %OutDat) {
	for my $item (sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$OutDat{$chr}}) {
		print O join("\t",$chr,@$item),"\n";
	}
}
close O;

__END__

perl ./getsyn.pl tigris2cat.out tigris2cat
perl ./getsyn.pl tigris2dog.out tigris2dog
perl ./getsyn.pl tigris2hum.out tigris2hum

join -1 4 -2 2 <(sort -k4 tigris2cat.plst) <(sort -k2 rec.pas) -o "1.1 1.2 1.3 1.4 2.3 2.7 2.8 2.9 2.10" -t $'\t' >plot2cat.lst
join -1 4 -2 2 <(sort -k4 tigris2dog.plst) <(sort -k2 rec.pas) -o "1.1 1.2 1.3 1.4 2.3 2.7 2.8 2.9 2.10" -t $'\t' >plot2dog.lst
join -1 4 -2 2 <(sort -k4 tigris2hum.plst) <(sort -k2 rec.pas) -o "1.1 1.2 1.3 1.4 2.3 2.7 2.8 2.9 2.10" -t $'\t' >plot2hum.lst

./dojoin.pl tigris2cat.plst 4 rec.npa 2 plot2cat
./dojoin.pl tigris2dog.plst 4 rec.npa 2 plot2dog
./dojoin.pl tigris2hum.plst 4 rec.npa 2 plot2hum

awk -F "\t" 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$10,$14,$15,$16,$17}' plot2cat.out |sort -k1.1 -k2.2n >plot2cat.nlst
awk -F "\t" 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$10,$14,$15,$16,$17}' plot2dog.out |sort -k1.1 -k2.2n >plot2dog.nlst
awk -F "\t" 'BEGIN {OFS="\t"} {print $1,$2,$3,$4,$10,$14,$15,$16,$17}' plot2hum.out |sort -k1.1 -k2.2n >plot2hum.nlst

perl -lane 'BEGIN {my %a} @x=split /\|/; $a{$x[-1]}=$_; END {print $a{$_} for sort keys %a}' plot2dog.chrorder > plot2dog.chrorder1
