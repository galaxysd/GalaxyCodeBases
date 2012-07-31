#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Data::Dump qw(ddx);

die "Usage: $0 <ec list> <out>\n" if @ARGV<1;
my $eclst=shift;
my $outfs=shift;

my $ReadLen = 101-5;
my $Eseq="CTGCAG";
my $EseqLen = length $Eseq;
my $EcutAt=5;

my $EfwdTerm5=1-$EcutAt+1;	# -3
my $EfwdTerm3=1-$EcutAt+$ReadLen;	# 92
my $ErevTerm3=0;	# 0
my $ErevTerm5=$ErevTerm3-$ReadLen+1;	# -95
my $Rfwd2Ec = -$EfwdTerm5;	# 3
my $Rrev2Ec = -$ErevTerm3;	# 0

# C|TGCAG
my $PosLeft = -$EfwdTerm3;	# -92
my $PosRight = -$ErevTerm5;	# 95;
my $PosECsft = $Rfwd2Ec;	# 3;
# - [x-$PosLeft,x+$PosECsft], Len=96
# + [x,x+$PosRight], Len=96
# A [x-$PosLeft,x+$PosRight], Len=95+92+1=186
my $minDis = $PosLeft + $PosRight;

open O,'>',$outfs or die "Error opening $outfs : $!\n";

my (%Markers,@ChrOrder);
open L,'<',$eclst or die;
while (<L>) {
	if (/^(#|$)/) {
		print O $_;
		next;
	}
	my ($chr,$pos,$strand,$mark,$count,$sampleCNT) = split /\t/;
	push @ChrOrder,$chr unless exists $Markers{$chr};
	push @{$Markers{$chr}{$pos}},[$strand,$mark,$count,$sampleCNT];
}
close L;

for my $chr (@ChrOrder) {
	my ($Fstrand,$Fmark,$Fcount,$FsampleCNT);
	for my $pos (keys %{$Markers{$chr}}) {
		next if @{$Markers{$chr}{$pos}} == 1;
		($Fstrand,$Fmark,$Fcount,$FsampleCNT)=(0,0,0,0);
		for (@{$Markers{$chr}{$pos}}) {
			my ($strand,$mark,$count,$sampleCNT)=@$_;
			if ($strand eq '+') {
				$Fstrand |= 1;
			} elsif ($strand eq '-') {
				$Fstrand |= 2;
			} else {
				$Fstrand = 3;
				$mark -= 10000;
			}
			$Fmark += $mark;
			$Fcount += $count;
			$FsampleCNT = $sampleCNT if $FsampleCNT < $sampleCNT;
		}
		if ($Fstrand == 1) {
			$Fstrand = '+';
		} elsif ($Fstrand == 2) {
			$Fstrand = '-';
		} else {
			$Fstrand = 'A';
			$Fmark += 10000;
		}
		$Markers{$chr}{$pos} = [[$Fstrand,$Fmark,$Fcount,$FsampleCNT]];
	}
}

my %tooNear;
for my $chr (@ChrOrder) {
	my @Pos = sort {$a <=> $b} keys %{$Markers{$chr}};
	my $lastPos=$Pos[0];
	for my $pos (@Pos) {
		print O join("\t",$chr,$pos,@{$Markers{$chr}{$pos}->[0]}),"\n";
		if ($pos != $lastPos and $pos-$lastPos < $minDis) {
			$tooNear{$chr}{$pos} += $Markers{$chr}{$pos}->[0][2];
		}
	}
}
close O;

ddx \%tooNear;

__END__
# radposbam_fix.pl:94: { scaffold413 => { 11969 => 1 }, scaffold574 => { 8197 => 349 } }
scaffold413	11967	A	10425.452997535	450	18
scaffold413	11969	-	0.801056595461641	1	1
scaffold574	8195	+	7	7	5
scaffold574	8197	A	10336.9968505454	349	17
