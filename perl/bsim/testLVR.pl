#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

use FindBin 1.51 qw($RealBin);
use lib "$RealBin";
use simVirusInserts;

sub InsertParts2RealPos($$$$$) {
	my ($r1fr,$PartsRef,$MappedChr,$simLMR,$VirSLR) = @_;
	my @Poses;
	for my $Part (@$PartsRef) {
		$Part =~ /^(\d+)(\w+)$/ or die "[$Part]";
		my ($slen,$type)=($1,$2);
		my $tmp;
		if ($type eq 'V') {
			if ( ($VirSLR->[0] eq '+' and $r1fr eq 'f') or ($VirSLR->[0] eq '-' and $r1fr eq 'r') ) {
				$tmp = $VirSLR->[1]; # LeftPos => not to: + $slen;
			} else {
				$tmp = 1+ $VirSLR->[2] - $slen;
			}
		} else {
			$tmp = 1+ $simLMR->[1] - $slen if $type eq 'L';	# `1+ $simLMR->[1] - $slen` can also be `1+ $simLMR->[0] + $innerPos` for 'L' && Read1。
			$tmp = 1+ $simLMR->[1] if $type eq 'R';	# LeftPos => not to add $slen
		}
		push @Poses,$tmp.$type."+$slen";
	}
	return @Poses;
}
sub getRealPos($$$$$$$$$$) {
	my ($r1fr,$innerPos, $InsertSize,$ReadLen,$VirSLR, $MappedChr,$rSMS,$r12,$strand,$simLMR)=@_;
	my @ret;	# 暂且忽略 $strand
	my $VirFrag = $VirSLR->[2] - $VirSLR->[1] +1;
	my ($FiveT,$ThreeT) = getInsertPos($r1fr,$innerPos,$InsertSize,$ReadLen,$r12);
	my @Parts = InsertPos2InsertParts($InsertSize,$ReadLen,$VirFrag, $FiveT,$ThreeT);
	my @Poses = InsertParts2RealPos($r1fr,\@Parts,$MappedChr,$simLMR,$VirSLR);
	ddx [$FiveT,$ThreeT,@Parts];
	return @Poses;
}

my $str = 'sf176_Ref_43898202_43898351_43898501_Vir_+_1878_1922_R_150_90';
$str =~ /^s([fr])(\d+)_Ref_(\d+)_(\d+)_(\d+)_Vir_([+-])_(\d+)_(\d+)_R_(\d+)_(\d+)$/ or die "$str";	# _([\d\|\-LVR]+)
my ($r1fr,$innerPos,$RefLeft,$RefMiddle,$RefRight,$VirStrand,$VirLeft,$VirRight,$InsertSize,$ReadLen) = ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10);

my ($r12R1,$r12R2,$r1SMS,$r2SMS)=(1,2);
my ($strandR1,$strandR2)=('+','+');
my ($refR1,$refR2)=('NA','NA');

my @Poses1 = getRealPos($r1fr,$innerPos,$InsertSize,$ReadLen, [$VirStrand,$VirLeft,$VirRight], $refR1,$r1SMS,$r12R1,$strandR1, [$RefLeft,$RefMiddle,$RefRight]);
my @Poses2 = getRealPos($r1fr,$innerPos,$InsertSize,$ReadLen, [$VirStrand,$VirLeft,$VirRight], $refR2,$r2SMS,$r12R2,$strandR2, [$RefLeft,$RefMiddle,$RefRight]);
