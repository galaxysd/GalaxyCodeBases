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

sub getInsertPos($$$$$) {
	my ($r1fr,$innerPos,$InsertSize,$ReadLen,$r12)=@_;
	my ($FiveT,$ThreeT)=(0,0);
	if (($r12 == 1 and $r1fr eq 'f') or ($r12 == 2 and $r1fr eq 'r')) {
		$FiveT = $innerPos;
		$ThreeT = $innerPos + $ReadLen;
	} elsif (($r12 == 2 and $r1fr eq 'f') or ($r12 == 1 and $r1fr eq 'r')) {
		my $InnerDis = $InsertSize - $ReadLen*2;
		$FiveT = $innerPos + $ReadLen + $InnerDis;
		$ThreeT = $innerPos + $InsertSize;
	} else {die 'Y';}
	return ($FiveT,$ThreeT);
}
sub InsertPos2PartLVR($$$$) {
	my ($r1fr,$Pos,$InsertSize,$VirFrag)=@_;
	my ($type,$lastingLen)=('',0);
	if ($r1fr eq 'f') {
		if ($Pos < $InsertSize) {
			$type = 'L';
			$lastingLen = $InsertSize - $Pos;
		} elsif ( $Pos < ($InsertSize + $VirFrag) ) {
			$type = 'V';
			$lastingLen = $VirFrag - ($Pos - $InsertSize);
		} else {
			$type = 'R';
			$lastingLen = 2*$InsertSize + $VirFrag - $Pos;
		}
	} elsif ($r1fr eq 'r') {
		if ( $Pos > ($InsertSize + $VirFrag) ) {
			$type = 'L';
			$lastingLen = $Pos - $InsertSize - $VirFrag;
		} elsif ( $Pos > $InsertSize ) {
			$type = 'V';
			$lastingLen = $Pos - $InsertSize;
		} else {
			$type = 'R';
			$lastingLen = $Pos;
		}
	}
	return ($type,$lastingLen);
}
sub Parts2List($$$$$$$$$) {
	my ($r1fr,$Pos,$type1,$lastingLen1,$type2,$lastingLen2,$ReadLen,$InsertSize,$VirFrag) = @_;
	my %Type2Int = (
		L => 1,
		V => 2,
		R => 3,
	);
	my %Int2Type = reverse %Type2Int;
	if ($type1 eq $type2) {
		return ( "${ReadLen}${type1}" );
	} elsif ( $Type2Int{$type1} < $Type2Int{$type2} ) {
		my @ret = ( "${lastingLen1}${type1}" );
		my $nextInt = $Type2Int{$type1};
		my $nextPos = $Pos;
		my $nextLL = $lastingLen1;
		my $nextType;
		my $RemainLen = $ReadLen;
		while ($nextInt < $Type2Int{$type2}) {
			++$nextInt;
			if ($r1fr eq 'f') {
				$nextPos += $nextLL;
			} else {
				$nextPos -= $nextLL;
			}
			$RemainLen -= $nextLL;
			($nextType,$nextLL) = InsertPos2PartLVR($r1fr,$nextPos,$InsertSize,$VirFrag);
			$nextLL=$RemainLen if $nextInt == $Type2Int{$type2};
			#ddx [$nextInt,$nextType,$RemainLen,$nextLL,1,$r1fr,$Pos,$type1,$lastingLen1,$type2,$lastingLen2,$ReadLen,$InsertSize,$VirFrag];
			push @ret,"${nextLL}$nextType";
		}
		return @ret;
	} else { die 'E'; }
}
sub InsertPos2InsertParts($$$$$$$$) {
	my ($InsertSize,$ReadLen,$VirLeft,$VirRight,$r1fr,$rSMS,$FiveT,$ThreeT)=@_;
	my $VirFrag = $VirRight - $VirLeft;
	my ($FiveSkip,$ThreeSkip)=(0,0);
	my ($type5,$lastingLen5) = InsertPos2PartLVR($r1fr,$FiveT,$InsertSize,$VirFrag);
	my ($type3,$lastingLen3) = InsertPos2PartLVR($r1fr,$ThreeT,$InsertSize,$VirFrag);
	my @Parts;
	if ($r1fr eq 'f') {
		$FiveSkip = $$rSMS[0];
		$ThreeSkip = $$rSMS[2];
		@Parts = Parts2List($r1fr,$FiveT,$type5,$lastingLen5,$type3,$lastingLen3,$ReadLen,$InsertSize,$VirFrag);
	} elsif ($r1fr eq 'r') {
		$FiveSkip = $$rSMS[2];
		$ThreeSkip = $$rSMS[0];
		@Parts = Parts2List($r1fr,$ThreeT,$type3,$lastingLen3,$type5,$lastingLen5,$ReadLen,$InsertSize,$VirFrag);
	}
	return @Parts;
}
sub InsertParts2RealPos() {
	;
}
sub getRealPos($$$$$$$$$$) {
	my ($r1fr,$innerPos, $InsertSize,$ReadLen,$VirLeft,$VirRight, $MappedChr,$rSMS,$r12,$strand)=@_;
	my ($FiveT,$ThreeT) = getInsertPos($r1fr,$innerPos,$InsertSize,$ReadLen,$r12);
	my @Parts = InsertPos2InsertParts($InsertSize,$ReadLen,$VirLeft,$VirRight,$r1fr,$rSMS, $FiveT,$ThreeT);
	#InsertParts2RealPos($MappedChr,$strand);
	return ($FiveT,$ThreeT);
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
		($R1Left,$R1Right) = getRealPos($r1fr,$innerPos,$InsertSize,$ReadLen,$VirLeft,$VirRight, $refR1,$r1SMS,$r12R1,$strandR1);
		($R2Left,$R2Right) = getRealPos($r1fr,$innerPos,$InsertSize,$ReadLen,$VirLeft,$VirRight, $refR2,$r2SMS,$r12R2,$strandR2);
		warn "$r1fr ($R1Left,$R1Right) $r12R1 ($R2Left,$R2Right) $r12R2 SamPos:$dat1[3],$dat2[3]\n";
	}
	close $bamfh;
	print STDERR ".\n";
}


ddx \%Stat;
__END__
./statVIsam.pl hs_ref_GRCh38.p2_chr18.mfa.gz HBV.AJ507799.2.fa.gz bam/*.sn.bam
