package main;
use strict;
use warnings;
use Digest::SHA;
use IPC::Open2;

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
		return '1CT';
	} elsif ($BaseCnt{'G'}<=$seqlen*$main::methly3BaseErrRate and $BaseCnt{'A'}>0) {
		return '2GA';
	} else {
		return '0Raw';
	}
}

# 由于include在main前面，所以可以在此统一计算打分矩阵。
my ($match,$methylmatch,$mismatch,$INDEL) = (2,1,-3,-4);
our (%MatrixO,%MatrixFct,%MatrixRga);
my @Bases = sort qw(A C G T);
# $ perl -e '@Bases = sort qw(A T C G); print join("|",@Bases),"\n"'
# A|C|G|T
for my $i (0 .. $#Bases) {
	for my $j ($i .. $#Bases) {
		my $str = $Bases[$i] . $Bases[$j];
		if ($i == $j) {
			$MatrixO{$str} = $MatrixFct{$str} = $MatrixRga{$str} = $match;
		} else {
			$MatrixO{$str} = $mismatch;
			if ($str =~ /^[CT]+$/) {
				$MatrixFct{$str} = $methylmatch;
			} else {
				$MatrixFct{$str} = $mismatch;
			}
			if ($str =~ /^[GA]+$/) {
				$MatrixRga{$str} = $methylmatch;
			} else {
				$MatrixRga{$str} = $mismatch;
			}
		}
	}
}
#ddx (\%MatrixO,\%MatrixFct,\%MatrixRga);
#	{ AA => 2, AC => -3, AG => -3, AT => -3, CC => 2, CG => -3, CT => -3, GG => 2, GT => -3, TT => 2 },
#	{ AA => 2, AC => -3, AG => -3, AT => -3, CC => 2, CG => -3, CT => 1, GG => 2, GT => -3, TT => 2 },
#	{ AA => 2, AC => -3, AG => 1, AT => -3, CC => 2, CG => -3, CT => -3, GG => 2, GT => -3, TT => 2 },

sub doAlign($$$) {
	my ($AssemHRef,$retHostARef,$retVirusARef) = @_;
	# 暂时只考虑第一条
	my $retHost = $$retHostARef[0];
	my $retVirus = $$retVirusARef[0];
	my @froDat = sort keys %{$AssemHRef};
	my $fro0 = shift @froDat;
	my $result = $AssemHRef->{$fro0}->[0];
	for my $fro (@froDat) {
		#$fro0 = mergeAln( $result,$AssemHRef->{$fro}->[0],$fro0,$fro );
		$result = mergeAln( $result,$AssemHRef->{$fro}->[0] );
	}
	#dynAln($retHost,$result);
	die;
}
sub dynAln($$$) {
	my ($ref,$query,$MatrixR) = @_;
	my @a=('',split //,$ref);
	my @b=('',split //,$query);
	my (@matrix,@path);
	my ($i,$j,$sc);
	$matrix[$_][0]=0 for (0 .. $#a);
	$matrix[0][$_]=0 for (0 .. $#b);
	$path[$#a][$#b] = 0;
	# dynamic recursion
	for ($i=1;$i<=$#a;$i++) {
		for ($j=1;$j<=$#b;$j++) {
			my $str = join('',sort ($a[$i],$b[$j]));
			$matrix[$i][$j] = $MatrixR->{$str} + $matrix[$i-1][$j-1];
			$path[$i][$j] = 1;	# case 0: i,j are aligned
			$sc=$matrix[$i-1][$j] + $INDEL;	# case 1: i aligned to -
			if ($sc > $matrix[$i][$j]) {
				$matrix[$i][$j] = $sc;
				$path[$i][$j] = 2;
			}
			$sc=$matrix[$i][$j-1] + $INDEL;	# case 2: j aligned to -
			if ($sc > $matrix[$i][$j]) {
				$matrix[$i][$j] = $sc;
				$path[$i][$j] = 3;
			}
			if ($matrix[$i][$j] < 0) {	# Local
				$matrix[$i][$j] = 0;
				$path[$i][$j] = 0;
			}
		}
	}
	# Traceback (with shadow matrix)
	my ($count,@ax,@ay,@as)=(0);	# first is 0.
	--$i;--$j;	# outside for, $i will bigger than $#a by 1 step
	my ($len,$mv,$mi,$mj)=(0);
	while ($i>=0 and $j>=0) {
		# Set Start Point
		#$mi=$i;$mj=$j;
		($mv,$mi,$mj)=@{&getmax(\@matrix,$i,$j)};	# Local
		# Traceback cycle for single str
		while ($mi>=0 and $mj>=0) {
			$as[$count] = $matrix[$mi][$mj] unless $as[$count];
			++$len;
			#&traceback($a[$mi],$b[$mj],$path[$mi][$mj],$ax[$count],$ay[$count]);
			my $path = $path[$mi][$mj];
			if (! defined $path) {	# $path == 0
				--$len;
				if ($len>0) {
					@{$ax[$count]} = reverse @{$ax[$count]};
					@{$ay[$count]} = reverse @{$ay[$count]};
					++$count;
					$len=0;
				}
				--$mi;--$mj;
				$i=$mi;$j=$mj;
				($mv,$mi,$mj)=@{&getmax(\@matrix,$i,$j)};	# Local
			} elsif ($path == 1) {
				push @{$ax[$count]},$a[$mi];push @{$ay[$count]},$b[$mj];
				--$mi;--$mj;
			} elsif ($path == 2) {
				push @{$ax[$count]},$a[$mi];push @{$ay[$count]},"-";
				--$mi;
			} elsif ($path == 3) {
				push @{$ax[$count]},"-";push @{$ay[$count]},$b[$mj];
				--$mj;
			} else {
				die;
			}
		}
	}
	for ($i=1;$i<=$count;$i++) {
		print "No. $i:\n";
		print 'X:[',@{$ax[$i-1]},"]\n",'Y:[',@{$ay[$i-1]},"]\n",'Score:',$as[$i-1],"\n";
	}
}
sub getmax($$$) {
	my ($arref,$ii,$ij)=@_;
	my ($mv,$mi,$mj,$i,$j)=(-1,-1,-1);
	for ($i=1;$i<=$ii;$i++) {
		for ($j=1;$j<=$ij;$j++) {
			if ($mv == $$arref[$i][$j]) {	# keep ($mi+$mj) max, to achieve max(mi^2+mj^2)
				if (($mi+$mj) < ($i+$j)) {
					#$mv=$$arref[$i][$j];
					$mi=$i;$mj=$j;
				}
			} elsif ($mv < $$arref[$i][$j]) {
				$mv=$$arref[$i][$j];
				$mi=$i;$mj=$j;
			}
		}
	}
	return [$mv,$mi,$mj];
}

sub mergeAln($$) {
	my ($ref,$query) = @_;
	my  $pid = open2( \*READER, \*WRITER, "./bin/alnmethly" );
	WRITER->autoflush(); # default here, actually
	my @Dat = ([$ref,undef],[$query,undef],[]);
	$Dat[2]->[0] = revcom($query);
	$_->[1] = guessMethyl($_->[0]) for @Dat;
	@Dat = sort { $a->[1] cmp $b->[1] } @Dat;
	print WRITER join("\n",map {$_->[0]} @Dat),"\n";
	while(<READER>) {
		chomp;
		if (/^Path(\d): ([IDMmR]+),(\d+)$/) {
			print "[$_] [$1] [$2] [$3]\n";
		}
	}
	print join("\n",map {$_->[0]} @Dat),"\n";
	return 1;
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
