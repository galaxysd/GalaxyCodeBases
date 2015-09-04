package main;
#use strict;
#use warnings;
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
	return 1 if exists $VirusChrIDs{$ChrA};
	return -1 if exists $VirusChrIDs{$ChrB};
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

sub mergeIn($$$$$) {
	my ($isHost,$rChrRange,$aRead,$rStore,$tid) = @_;
	my ($reflen,$readlen) = cigar2poses($aRead->[5]);
	my $thisehPos = $aRead->[3]+$reflen;
	my $ret;	# 1 -> again, 0 -> merged
	if (exists $rChrRange->{$aRead->[2]}) {
		if ($aRead->[3] <= $rChrRange->{$aRead->[2]}->[1]) {
			#die unless $thisehPos >= $rChrRange->{$aRead->[2]}->[1];
			$rChrRange->{$aRead->[2]}->[1] = $thisehPos if $thisehPos >= $rChrRange->{$aRead->[2]}->[1];
			$ret = 0;
		} else {
			if ($isHost) {
				$ret = 1;
			} else {
				$ret = 0;
				# assume overlap
				if ($aRead->[3] <= $rChrRange->{$aRead->[2]}->[0]) {
					$rChrRange->{$aRead->[2]}->[0] = $aRead->[3];
				}
				if ($thisehPos >= $rChrRange->{$aRead->[2]}->[1]) {
					$rChrRange->{$aRead->[2]}->[1] = $thisehPos;
				}
			}
		}
	} else {
		if (keys %{$rChrRange}) {
			if ($isHost) {
				$ret = 1;
			} else {
				$ret = 0;
				$rChrRange->{$aRead->[2]} = [ $aRead->[3],$thisehPos ];
			}
		} else {
			$rChrRange->{$aRead->[2]} = [ $aRead->[3],$thisehPos ];
			$ret = 0;
		}
	}
	unless ($ret) {
		push @{$rStore},$tid;
	}
	return $ret;
}
sub formatChrRange($) {
	my ($rChrRange) = @_;
	my @ret;
	for my $c (sort keys %{$rChrRange}) {
		push @ret, "${c}:".$rChrRange->{$c}->[0].'-'.$rChrRange->{$c}->[1];
	}
	return join(',',@ret);
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
