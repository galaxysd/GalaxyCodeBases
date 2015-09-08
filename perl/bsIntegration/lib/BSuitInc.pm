package main;
use strict;
use warnings;
use Digest::SHA;

sub getFilesHash(@) {
	my $fileStr = join(',',@_);
	return Digest::SHA::sha1_base64($fileStr);
}

sub getRef2char($$) {
	my ($HostRefName,$VirusRefName)=@_;
	my $HostChar = substr $HostRefName,0,1;
	my $VirusChar = substr $VirusRefName,0,1;
	my $tlen = length($VirusRefName) -1;
	my $i=0;
	while ( ($VirusChar eq $HostChar) and (++$i<=$tlen) ) {
		my $tmpChar = substr $VirusRefName,$i,1;
		$VirusChar = $tmpChar if $tmpChar=~/\w/;
		warn "$i - $tmpChar $VirusChar -\n";
	}
	$VirusChar = 'V' if $VirusChar eq $HostChar;
	return $HostChar.$VirusChar;
}

sub sortChrPos($$) {
	my ($ChrA,$PosA) = split /\t/,$_[0];
	my ($ChrB,$PosB) = split /\t/,$_[1];
	if ($ChrA eq $ChrB) {
		return $PosA <=> $PosB;
	}
	return 1 if exists $main::VirusChrIDs{$ChrA};
	return -1 if exists $main::VirusChrIDs{$ChrB};
	$ChrA cmp $ChrB ||
	$PosA <=> $PosB;
}

sub cigar2poses($) {
	my ($cigar) = @_;
	my @cigar = $cigar =~ /(\d+)(\w)/g;
	my ($reflen,$maxM,$readlen)=(0,0,0);
	while (@cigar) {
		my ($len,$op) = splice(@cigar,0,2);
		if ($op eq 'M') {
			$reflen += $len;
			$readlen += $len;
			$maxM = $len if $maxM < $len;
			next;
		}
		$reflen += $len if $op eq 'D';
		$readlen += $len if $op eq 'I';
	}
	return ($reflen,$readlen);
}

sub mergeIn($$$$$$) {
	my ($isHost,$rChrRange,$aRead,$rStore,$cid,$i) = @_;
	my ($reflen,$readlen) = cigar2poses($aRead->[5]);
	my $thisehPos = $aRead->[3]+$reflen;
	my $ret;	# 1 -> again, 0 -> merged
	if ($i > 1 or $isHost == 0 or scalar(keys %{$rChrRange})==0) {
		# same PE or assume overlap on Virus
		$ret = 0;
		if (keys %{$rChrRange} and exists $rChrRange->{$aRead->[2]}) {
			if ($aRead->[3] <= $rChrRange->{$aRead->[2]}->[0]) {
				$rChrRange->{$aRead->[2]}->[0] = $aRead->[3];
			}
			if ($thisehPos >= $rChrRange->{$aRead->[2]}->[1]) {
				$rChrRange->{$aRead->[2]}->[1] = $thisehPos;
			}
		} else {
			$rChrRange->{$aRead->[2]} = [ $aRead->[3],$thisehPos,0 ];
		}
	} elsif (exists $rChrRange->{$aRead->[2]}) {
		if ($aRead->[3] <= $rChrRange->{$aRead->[2]}->[1]) {	# overlap
			#die unless $thisehPos >= $rChrRange->{$aRead->[2]}->[1];
			$rChrRange->{$aRead->[2]}->[1] = $thisehPos if $thisehPos >= $rChrRange->{$aRead->[2]}->[1];
			$ret = 0;
		} else {
			$ret = 1;
		}
	} else {	# different HostChr
		$ret = 1;
	}
	unless ($ret) {
		my $tid = join("\n",$cid,$i);
		push @{$rStore},$tid;
		++$rChrRange->{$aRead->[2]}->[2];
	}
	return $ret;
}
sub formatChrRange($) {
	my ($rChrRange) = @_;
	my @ret;
	for my $c (sort keys %{$rChrRange}) {
		push @ret, "${c}:".$rChrRange->{$c}->[0].'-'.$rChrRange->{$c}->[1].':'.$rChrRange->{$c}->[2];
	}
	return join(',',@ret);
}

sub revcom($) {
	my $str = $_[0];
	$str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
	my $rev = reverse $str;
	$rev =~ tr/[](){}<>/][)(}{></;
	return $rev;
}

sub guessMethyl($) {
	my ($seq) = @_;
	my %BaseCnt=(
		A => 0, G => 0, C => 0, T => 0,
		N => 0,
	);
	for (split //,$seq) {
		++$BaseCnt{uc $_};
	}
	my $seqlen = length($seq) - $BaseCnt{'N'};
	return 'N' if $seqlen == 0;
	my @Cnts = sort { $BaseCnt{uc $b} <=> $BaseCnt{uc $a} } keys %BaseCnt;
	#ddx [\@Cnts,\%BaseCnt];
	if ($BaseCnt{'C'}<=$seqlen*$main::methly3BaseErrRate and $BaseCnt{'T'}>0) {
		return 'CT';
	} elsif ($BaseCnt{'G'}<=$seqlen*$main::methly3BaseErrRate and $BaseCnt{'A'}>0) {
		return 'GA';
	} else {
		return 'Raw';
	}
}

sub warnFileExist(@) {
	my %NotFound;
	for (@_) {
		++$NotFound{$_} unless -f $_;
	}
	my @NF = sort keys %NotFound;
	if (@NF > 0) {
		warn "[!!!] File NOT Found:[",join('],[',@NF),"]\n";
	}
	#warn "[Debug] @_\n";
	return join(' ',@_);
}
# http://perldoc.perl.org/perlfaq4.html#How-do-I-expand-function-calls-in-a-string%3f

1;
