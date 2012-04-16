#!/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120412
=cut
use strict;
use warnings;
use Data::Dump qw(ddx);

die "Usage: $0 <in file> <out file>\n" if @ARGV<2;
my $in=shift;
my $out=shift;

my @Indi;
my ($indi,$site)=(0,0);

open I,'<',$in or die;
while (<I>) {
	next if /^#/;
	chomp;
	my (undef,undef,$id,@GT)=split /\t/;
	$Indi[$id] = \@GT;
	die if $indi && $indi != @GT;
	$indi = scalar @GT;
}
$site = $#Indi;
push @{$Indi[0]},0 for (0 .. $indi-1);
print "Loaded: $site SNP sites * $indi Sperms.\n";

my @a=split //,('0' x $site);
my @b=split //,('1' x $site);
my (@preStats,$preStatsAt);
push @preStats,[split //,('0' x $site)] for (1..($site-1));	# Write here for calloc in C.
$preStatsAt=0;
my (@matrix,@path);	# 2D array in perl can be malloced auto.ly
my ($i,$j,$sc,$mstep);
my %COLOR=('1'=>'32','2'=>'33','3'=>'36','0'=>'0');

sub scoring($$$) {	# Recombin_Count >> 0, so 0 is reserved for NULL.
	my ($i,$LR,$lastScore)=@_;
	die if ($mstep-$preStatsAt > 1);	# @preStats must be the latest one. Debug only.
	my $RC=0;
	for my $n (1..$indi) {
		if ( $Indi[$mstep][$n] != ($LR^$preStats[$i][$n]) ) {	# remember the higher precedence of '!=' to '^'
			++$RC;
		}
	}
	return $lastScore+$RC;
}
sub updatePreStats($$) {
	my ($i,$LR)=@_;
	die if ($mstep-$preStatsAt > 1);	# @preStats must be the latest one. Debug only.
	if ($LR==3) {	# remain old status on half/half.
		return;
	} elsif ($LR == 1) {
		$LR = 0;
	} elsif ($LR == 2) {
		$LR = 1;
	} else {die "Impossible for \$LR>3.";}
	for my $n (1..$indi) {
		$preStats[$i][$n] = ($Indi[$mstep][$n]^$LR);
	}
	return;
}

# Initialization
$matrix[0][0]=0;
$path[0][0]=0;
for ($i=0;$i<=$site;$i++) {
	for ($j=0;$j<=$site;$j++) {
		$matrix[$i][$j] = 0;
	}
}
#for ($j=1;$j<=$#b;$j++) {$matrix[0][$j] = $j*$INDEL;$path[0][$j]=3;}
=pod
 ij                                    j => R
         00, 0, -      01, 1, R      02, 2, RR     03, 3, RRR    04, 4, RRRR
         10, 1, L      11, 2, LR/RL  12, 3, NNN    13, 4, NNNN   (14, 5)
 i => L  20, 2, LL     21, 3, NNN    22, 4, NNNN   (23, 5)       (24, 6)
         30, 3, LLL    31, 4, NNNN   (32, 5)       (33, 6)       (34, 7)
         40, 4, LLLL   (41, 5)       (42, 6)       (43, 7)       (44, 8)
=cut
for ($mstep=1;$mstep<=$site;$mstep++) {	# starts from 01 & 10.
	#updatePreStats();
	for ($i=0;$i<=$mstep;$i++) {
		$j = $mstep-$i;	# only half matrix needed.
		my ($sc1,$sc2)=(0,0);
		if ($j>0) {	# has Up, can be L
			$sc1=scoring($i,0,$matrix[$i][$j-1]);
		}
		if ($i>0) {	# has Left, can be R
			$sc2=scoring($i,1,$matrix[$i-1][$j]);
		}
		if ($sc1) {
			if ($sc2) {	# both $sc1 and $sc2
				if ($sc1 < $sc2) {
					$path[$i][$j]=1;
					$sc = $sc1;
				} elsif ($sc1 == $sc2) {
					$path[$i][$j]=3;	# half L => 0.5, no difference.
					$sc = $sc1;
				} else {	# $sc1 > $sc2
					$path[$i][$j]=2;
					$sc = $sc2;
				}
			} else {	# only $sc1
				$path[$i][$j]=1;
				$sc = $sc1;
			}
		} elsif ($sc2) {	# only $sc2
			$path[$i][$j]=2;
			$sc = $sc2;
		} else {
			die "Impossible for $i==$j==0.";
			$path[$i][$j]=0;
			$sc = 0;
		}
		$matrix[$i][$j] = $sc;
		updatePreStats($i,$path[$i][$j]);
	}
	++$preStatsAt;
}

for ($mstep=1;$mstep<=$site;$mstep++) {	# starts from 01 & 10.
	for ($i=0;$i<=$mstep;$i++) {
		$j = $mstep-$i;
		print $matrix[$i][$j],',',$path[$i][$j],"\t";
	}
	print "\n";
}
