#!/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120412
=cut
use strict;
use warnings;
#use Data::Dump qw(ddx);

die "Usage: $0 <in file> <out file>\n" if @ARGV<2;
my $in=shift;
my $out=shift;

my (@Indi,@Parents);
my ($indi,$site)=(0,0);

open I,'<',$in or die;
while (<I>) {
	next if /^#/;
	chomp;
	my ($parent,undef,$id,@GT)=split /\t/;
	$Indi[$id] = \@GT;	# @Indi[SNP][0..$indi-1]
	push @Parents,$parent;
	die if $indi && $indi != @GT;
	$indi = scalar @GT;
}
$site = $#Indi;	# count of data lines, 1st is 0. @Indi[1..SNP][0..$indi-1]
push @{$Indi[0]},0 for (0 .. $indi-1);	# also set @Indi[SNP==0][0..$indi-1]=0
print "Loaded: $site SNP sites * $indi Sperms.\n";

my @a=split //,('0' x $site);
my @b=split //,('1' x $site);
my (@matrix,@path);	# 2D array in perl can be malloced auto.ly. @matrix[0..SNP][0..SNP]
my ($i,$j,$sc,$mstep);
my %COLOR=('1'=>'32','2'=>'33','3'=>'36','0'=>'0');

sub scoring($$$) {	# Recombin_Count >> 0, so 0 is reserved for NULL.
	my ($i,$LR,$lastGrid)=@_;
	my ($lasti,$lastj)=@$lastGrid;
	my $lastLR;
	if ($path[$lasti][$lastj] == 1) {
		$lastLR = 0;
	} else { $lastLR = 1; }
	my $RC=0;
	for my $n (0..$indi-1) {
		if ( ($LR^$Indi[$mstep][$n]) != ($lastLR^$Indi[$mstep-1][$n]) ) {	# remember the higher precedence of '!=' to '^'
			++$RC;
		}
	}
	return $matrix[$lasti][$lastj]+$RC;
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
	for ($i=0;$i<=$mstep;$i++) {
		$j = $mstep-$i;	# only half matrix needed.
		my ($sc1,$sc2)=(0,0);	# better to set NULL as -1.
		if ($j>0) {	# has Up, can be L
			$sc1=scoring($i,0,[$i,$j-1]);
		}
		if ($i>0) {	# has Left, can be R
			$sc2=scoring($i,1,[$i-1,$j]);
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
		} else {	# $sc1==$sc2==0
			#die "Impossible for $i==$j==0. $sc1,$sc2.";
			$path[$i][$j]=0;	# for path, 00 <=> 11
			$sc = 0;
		}
		$matrix[$i][$j] = $sc;
	}
}
my ($min,$mini,$minj)=($matrix[0][$site],0,$site);
for ($i=1;$i<=$site;$i++) {	# init value is from $i==0
	$j = $site-$i;
	if ($min > $matrix[$i][$j]) {
		$min=$matrix[$i][$j];
		($mini,$minj)=($i,$j);
	}
}
# Traceback
my ($sumxorP,@resultP,@xorP)=(0);
while ($mini>=0 and $minj>=0 and $mini+$minj) {
	push @resultP,$path[$mini][$minj];
	if ($path[$mini][$minj] & 1) {
		--$minj;
	} else { --$mini; }
}
@resultP = reverse @resultP;
for (0..$#resultP) {
	my $t = ($resultP[$_]^$Parents[$_])&1;
	push @xorP,$t;
	$sumxorP += $t;
}

print "Parents: ",join(',',@Parents),"\n";
print "Results: ",join(',',@resultP),"\n";
print "bitXORs: ",join(',',@xorP),"\nSumbXOR: $sumxorP\n";
print '-' x 78,"\nColor:";
for (sort keys %COLOR) {
	print "\033[",$COLOR{$_},";1m $_";
}
print "\033[0;0m\n";

open O,'>',$out or die;
print O "Loaded: $site SNP sites * $indi Sperms.\n";

for ($i=0;$i<=$site;$i++) {
	for ($j=0;$j<=$site-$i;$j++) {
		print "\033[",$COLOR{$path[$i][$j]},";1m";
		printf("%4d ",$matrix[$i][$j]);
		print O $matrix[$i][$j],',',$path[$i][$j],"\t";
	}
	print "\033[0;0m\n";
	print O "\n";
}

print O '-' x 78,"\n";
print O "Parents: ",join(',',@Parents),"\n";
print O "Results: ",join(',',@resultP),"\n";
print O "bitXORs: ",join(',',@xorP),"\nSumbXOR: $sumxorP\n";
close O;
