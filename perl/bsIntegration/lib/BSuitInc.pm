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
