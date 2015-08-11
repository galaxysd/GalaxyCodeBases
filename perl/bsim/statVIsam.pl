#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

use FindBin 1.51 qw($RealBin);
use lib "$RealBin";
use simVirusInserts;

die "Usage: $0 <Host> <Virus> <bam sort-n> [more bam files]\n" if @ARGV <3;

my $Reff = shift;
my $Virf = shift;

#$Reff='hs_ref_GRCh38.p2_chr18.mfa.gz';
#$Virf='HBV.AJ507799.2.fa.gz';

my $Refh = openfile($Reff);
my $RefID = getRefChr1stID($Refh);
close $Refh;
my $Virfh = openfile($Virf);
my $VirID = getRefChr1stID($Virfh);
close $Virfh;
warn "[!]Ref:[$RefID], Virus:[$VirID].\n";

sub InsertParts2RealPos($$$$$$) {
	my ($r1fr,$rSMS,$PartsRef,$MappedChr,$simLMR,$VirSLR) = @_;
	my ($FiveSkip,$ThreeSkip)=(0,0);
	#ddx $PartsRef;
	if ($r1fr eq 'f') {
		$FiveSkip = $$rSMS[0];
		$ThreeSkip = $$rSMS[2];
	} elsif ($r1fr eq 'r') {
		$FiveSkip = $$rSMS[2];
		$ThreeSkip = $$rSMS[0];
	}
	my @Poses;
	for my $Part (@$PartsRef) {
		$Part =~ /^(\d+)(\w+)$/ or die "[$Part]";
		my ($slen,$type)=($1,$2);
		my $tmp;
		if ($type eq 'V') {
			if ( ($VirSLR->[0] eq '+' and $r1fr eq 'f') or ($VirSLR->[0] eq '-' and $r1fr eq 'r') ) {
				$tmp = 1+ $VirSLR->[1]; # LeftPos => not to: + $slen;
			} else {
				$tmp = $VirSLR->[2] - $slen;
			}
		} else {
			$tmp = 1+ $simLMR->[1] - $slen if $type eq 'L';
			$tmp = 1+ $simLMR->[1] if $type eq 'R';	# LeftPos => not to add $slen
		}
		push @Poses,$tmp.$type."+$slen";
	}
	return @Poses;
}
sub getRealPos($$$$$$$$$$) {
	my ($r1fr,$innerPos, $InsertSize,$ReadLen,$VirSLR, $MappedChr,$rSMS,$r12,$strand,$simLMR)=@_;
	my @ret;	# 暂且忽略 $strand
	my $VirFrag = $VirSLR->[2] - $VirSLR->[1];
	my ($FiveT,$ThreeT) = getInsertPos($r1fr,$innerPos,$InsertSize,$ReadLen,$r12);
	my @Parts = InsertPos2InsertParts($InsertSize,$ReadLen,$VirFrag,$r1fr, $FiveT,$ThreeT);
	my @Poses = InsertParts2RealPos($r1fr,$rSMS,\@Parts,$MappedChr,$simLMR,$VirSLR);
	ddx [$FiveT,$ThreeT,@Parts];
	return @Poses;
}

my %Stat;
for my $bamin (@ARGV) {
	print STDERR "[!]Reading [$bamin] ";
	my $bamfh = openfile($bamin);
	while (<$bamfh>) {
		my @dat1 = split /\t/;
		$_=<$bamfh>;
		my @dat2 = split /\t/;
		die "[x]Read1 & Read2 not match ! [$dat1[0]] ne [$dat2[0]]\n" if $dat1[0] ne $dat2[0];
		# sf0_Ref_2707868_2708068_2708268_Vir_-_5629_5731
		$dat1[0] =~ /^s([fr])(\d+)_Ref_(\d+)_(\d+)_(\d+)_Vir_([+-])_(\d+)_(\d+)_R_(\d+)_(\d+)$/ or die "$dat1[0]";
		my ($r1fr,$innerPos,$RefLeft,$RefMiddle,$RefRight,$VirStrand,$VirLeft,$VirRight,$InsertSize,$ReadLen) = ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10);
		#print "$dat1[0] $innerPos,$RefLeft,$RefMiddle,$RefRight,$VirStrand,$VirLeft,$VirRight\n";
		next if $r1fr eq 'r';	# 封印反向Reads。
		my $r1SMS = cigar2SMS($dat1[5]);
		my $r2SMS = cigar2SMS($dat2[5]);
		print "$dat1[0] [",join('.',@$r1SMS),"][",join('.',@$r1SMS),"]\n";
		my ($r12R1,$r12R2)=(0,0);
		if (($dat1[1] & 0x40) and ($dat2[1] & 0x80) ) {
			$r12R1 = 1; $r12R2 = 2;
		} elsif (($dat1[1] & 0x80) and ($dat2[1] & 0x40) ) {
			$r12R1 = 2; $r12R2 = 1;
		} else {die 'X';}
		my ($strandR1,$strandR2)=('+','+');
		$strandR1 = '-' if $dat1[1] & 0x10;
		$strandR2 = '-' if $dat2[1] & 0x10;
		my ($refR1,$refR2)=('NA','NA');
		if ($dat1[2] eq $RefID) {
			$refR1 = 'Host';
		} elsif ($dat1[2] eq $VirID) {
			$refR1 = 'Virus';
		} else {$refR1 = "Other:$dat1[2]";}
		if ($dat2[2] eq $RefID) {
			$refR2 = 'Host';
		} elsif ($dat2[2] eq $VirID) {
			$refR2 = 'Virus';
		} else {$refR2 = "Other:$dat2[2]";}
		my ($R1Left,$R1Right,$R2Left,$R2Right)=(0,0,0,0);
		my @Poses1 = getRealPos($r1fr,$innerPos,$InsertSize,$ReadLen, [$VirStrand,$VirLeft,$VirRight], $refR1,$r1SMS,$r12R1,$strandR1, [$RefLeft,$RefMiddle,$RefRight]);
		my @Poses2 = getRealPos($r1fr,$innerPos,$InsertSize,$ReadLen, [$VirStrand,$VirLeft,$VirRight], $refR2,$r2SMS,$r12R2,$strandR2, [$RefLeft,$RefMiddle,$RefRight]);
		print "(@Poses1)->@dat1[3,5] $refR1\n(@Poses2)->@dat2[3,5] $refR2\n";
	}
	close $bamfh;
	print STDERR ".\n";
}


ddx \%Stat;
__END__
./statVIsam.pl chr18.fa HBV.AJ507799.2.fa.gz bam/Fsimout_m7C.sn.bam
bam/*.sn.bam

输出：
sf94_Ref_30669853_30670003_30670153_Vir_+_88313_88358_R_150_90 [0.56.34][0.56.34]
(30669948L+56 88348V+34)->30669948 56M34S Host
(88355V+41 30670053R+49)->30670004 41S49M Host
括号内是模拟的位置，箭头右边是bam文件的内容。
第一条是Read1，第二条是Read2.

30669948L+56：模拟在左侧，人的30669948处，然后有56M
88348V+34 后面跟着的病毒序列是88348开始的，有34M
可以看出，bwa把它比对到长的地方了。

负链的上周折腾了好几天，你找人帮忙写下吧。现在的结果是右端点的坐标。
