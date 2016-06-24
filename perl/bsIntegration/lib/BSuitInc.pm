package main;
use strict;
use warnings;
use Digest::SHA;
use IPC::Open2;
use Galaxy::Data;

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
	#ddx $AssemHRef;
	my $fro0 = shift @froDat;
	my $result = [$AssemHRef->{$fro0}->[0]->[1]];
	for my $fro (@froDat) {
		#$fro0 = mergeAln( $result,$AssemHRef->{$fro}->[0],$fro0,$fro );
		$result = mergeAln( $result->[0],$AssemHRef->{$fro}->[0]->[1] );
	}
	#ddx $retHost,$result;
	my $HostResult = mergeAln( $retHost->[1],$result->[0] );
	#ddx $HostResult;
	my @RefInfo = split /[:-]/,$retHost->[0];
	my ($RefCut,$Left,$Insert,@Rcuts,@Vcuts)=(0);
	if ($HostResult->[1] =~ /^(D*)([MmR]+)(I+)/) {
		(undef,$Left,$Insert) = ($1,$2,$3);
		push @Rcuts,($RefInfo[1] + length($Left)+1);
	}
	if ($HostResult->[1] =~ /^(I*)([MmR]+)(D+)/) {
		($Left,$Insert,undef) = ($1,$2,$3);
		push @Rcuts,($RefInfo[1] + length($Left) + length($Insert)+1);
	}
	$RefCut = join(';',@Rcuts);
	my $VirusResult = mergeAln( $retVirus->[1],$result->[0] );
	my @VirusInfo = split /[:-]/,$retVirus->[0];
	my ($VirCut,$VirLeft,$VirInsert)=(0,0,0);
	if ($VirusResult->[1] =~ /^(D*)([MmR]+)(I+)/) {
		(undef,$VirLeft,$VirInsert) = ($1,$2,$3);
		push @Vcuts,($VirusInfo[1] + length($VirLeft)+1);
	}
	if ($VirusResult->[1] =~ /^(I*)([MmR]+)(D+)/) {
		($Left,$Insert,undef) = ($1,$2,$3);
		push @Vcuts,($VirusInfo[1] + length($Left) + length($Insert)+1);
	}
	#ddx $VirusResult;
	$VirCut = join(';',@Vcuts);
	return [$RefInfo[0],$RefCut,$VirusInfo[0],$VirCut,length($VirInsert)];
}
sub dynAln($$$) { # 废弃 {
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
sub getmax($$$) { # 废弃 {
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
	my  $pid = open2( \*READER, \*WRITER, "$RealBin/bin/alnmethly" );
	WRITER->autoflush(); # default here, actually
	my @Dat = ([$query,undef],[revcom($query),undef]);
	$_->[1] = guessMethyl($_->[0]) for @Dat;
	@Dat = sort { $a->[1] cmp $b->[1] } @Dat;
	unshift @Dat,[$ref,guessMethyl($ref)];
	for (@Dat) {
		if ($_->[1] eq '1CT') {
			$_->[0] =~ s/[CT]/Y/ig;
		} elsif ($_->[1] eq '2GA') {
			$_->[0] =~ s/[GA]/R/ig;
		}
		$_->[0]='N' if $_->[1] eq 'N';	# 'N' for empty seq.
	}
	#ddx \@Dat;
	print WRITER join("\n",map {$_->[0]} @Dat),"\n";
	my %Result;
	while(<READER>) {
		chomp;
		if (/^Path(\d): ([IDMmR]+),(\d+)$/) {
			#print "[$_] [$1] [$2] [$3]\n";
			$Result{$1} = [$2,$3];
		}
	}
	waitpid( $pid, 0 );
	close READER;
	my @Resu = sort { $Result{$b}->[1] <=> $Result{$a}->[1] } keys %Result;
	#print join("\n",map {$_->[0]} @Dat),"\n";
	#ddx $Result{$Resu[0]},$Resu[0];
	my @ResDat = split //,$Result{$Resu[0]}->[0];
	my @QuaryDat = split //,$Dat[$Resu[0]]->[0];
	my @Refdat = split //,$ref;
	#print "R:@Refdat\nQ:@QuaryDat\nA:@ResDat\n";
	my @AlnDat;
	my ($Refp,$Querp) = (0,0);
	for my $p (0 .. $#ResDat) {
		if ($ResDat[$p] =~ /^[DMmR]$/) {
			++$Querp if $ResDat[$p] eq 'D';
			my $theBases = $IUB{$Refdat[$p-$Refp]} or die;
			push @$theBases,@$theBases if $#$theBases==0;
			++$AlnDat[$p]->{$_} for @$theBases;
		}
		if ($ResDat[$p] =~ /^[MmR]$/) {
			my $theBases = $IUB{$QuaryDat[$p-$Querp]};
			push @$theBases,@$theBases if $#$theBases==0;
			++$AlnDat[$p]->{$_} for @$theBases;
		} elsif ($ResDat[$p] eq 'I') {
			++$Refp;
			my $theBases = $IUB{$QuaryDat[$p-$Querp]};
			push @$theBases,@$theBases if $#$theBases==0;
			++$AlnDat[$p]->{$_} for @$theBases;
		}
	}
	#ddx \@AlnDat;
	my $reseq;
	for (@AlnDat) {
		my %t = %{$_};
		my @k = sort { $t{$a} <=> $t{$b} } keys %t;
		my $iub;
		if (@k>1 and $t{$k[0]} == $t{$k[1]}) {
			$iub = join('',sort @k[0,1]);
		} else {
			$iub = $k[0];
		}
		$reseq .= $REV_IUB{$iub};
	}
	return [$reseq,$Result{$Resu[0]}->[0]];
}

sub warnFileExist(@) {
	my %NotFound;
	for (@_) {
		next if /^\s*$/;
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

=pod
#define BAM_CIGAR_STR   "MIDNSHP=XB"
/* bam_cigar_type returns a bit flag with:
 *   bit 1 set if the cigar operation consumes the query
 *   bit 2 set if the cigar operation consumes the reference
 *
 * For reference, the unobfuscated truth table for this function is:
 * BAM_CIGAR_TYPE  QUERY  REFERENCE
 * --------------------------------
 * BAM_CMATCH      1      1
 * BAM_CINS        1      0
 * BAM_CDEL        0      1
 * BAM_CREF_SKIP   0      1
 * BAM_CSOFT_CLIP  1      0
 * BAM_CHARD_CLIP  0      0
 * BAM_CPAD        0      0
 * BAM_CEQUAL      1      1
 * BAM_CDIFF       1      1
 * BAM_CBACK       0      0
 * --------------------------------
 */
=cut
sub bam_cigar2rlen($) {
	my @cigar = $_[0] =~ /(\d+\w)/g;
	my $pos = 0;
	for (@cigar) {
		if (/(\d+)([MDN=X])/) {
			$pos += $1;
		}
	}
	return $pos;
}
sub bam_cigar2qlen($) {
	my @cigar = $_[0] =~ /(\d+\w)/g;
	#ddx \@cigar;
	my $pos = 0;
	for (@cigar) {
		#/(\d+)([MIS=X])/ or die $_[0],",$_";
		if (/(\d+)([MIS=X])/) {
			$pos += $1;
		}
	}
	return $pos;
}
sub getSeqCIGAR($$) {
	my ($minLeft,$i) = @_;
	my $firstSC=0;
	my $deltraPos = $i->[3] - $minLeft;
	my @cigar = $i->[5] =~ /(\d+[MIDNSHP=XB])/g;
	my $pos = $i->[3];
	my $offset = $pos - $minLeft;
	#my $cigar2qlen = bam_cigar2qlen($i->[5]);
	my ($cursorQ,$cursorR) = (0,0);
	my $seqCIGAR = 'B' x $offset;
	#die $offset;
	for (@cigar) {
		/(\d+)([MIDNSHP=XB])/ or die;
		if ($2 eq 'S') {
			if ($cursorQ == 0) {
				substr $seqCIGAR,-$1,$1,$2 x $1;
				$firstSC = $1;
			} else {
				$cursorQ += $1;
				$seqCIGAR .= $2 x $1;
			}
		} elsif ($2 eq 'M') {
			$cursorQ += $1;
			$seqCIGAR .= $2 x $1;
			for my $p ( ($cursorQ-$1) .. ($cursorQ-1) ) {
				my $refpos = $pos + $cursorR + $p;
				#$ReadsMPos2Ref{$p} = $refpos;
			}
			$cursorR += $1;
		} elsif ($2 eq 'I') {
			$cursorQ += $1;
			$seqCIGAR .= $2 x $1;
		} elsif ($2 eq 'D') {
			$cursorR += $1;
		} else {die;}
	}
	return ($firstSC,$seqCIGAR);
}
sub grepmerge($) {
my $DEBGUHERE = 0;
	my ($minLeft,$maxLS,@clipReads,%relPoses,%relPosesFR,$chr)=(5000000000,0);
	for my $i (@{$_[0]}) {
		my @cigar = $i->[5] =~ /(\d+[MIDNSHP=XB])/g;
		my $flag = 0;
		#ddx \@cigar; die $i->[5];
		for (@cigar) {
			if (/(\d+)S/) {
				$flag = 1 if $1 >= $main::minSoftClip;
			}
		}
		if ($flag) {
			unless (defined $chr) {
				$chr = $i->[2];
			} else {
				next if $chr ne $i->[2];
			}
#print " $i->[2]:$i->[3]:$i->[-1]\n" if $DEBGUHERE;
			push @clipReads,$i;
			my $left = $i->[3];
			if ($cigar[0] =~ /(\d+)S/) {
				$maxLS = $1 if $maxLS < $1;
				$left -= $1;
			}
			$minLeft = $left if $minLeft > $left;
		}
	}
	return ([],[]) unless @clipReads;
	#ddx \@clipReads;
	#my (@b2pCIGAR,@b2pDeriv);
	my @ReadsCIGAR;
print "$minLeft\n" if $DEBGUHERE;
	for my $i (@clipReads) {
		#my %ReadsMPos2Ref;
		my ($firstSC,$seqCIGAR) = getSeqCIGAR($minLeft,$i);
		my $offset = $i->[3] - $minLeft;
		my @Poses = getDeriv("M","B",$seqCIGAR);
		for (@Poses) {
			++$relPosesFR{$_};
			if ($_ >= 0) {
				++$relPoses{$_+0.5};
			} else {
				++$relPoses{-$_-0.5};
			}
		}
if ($DEBGUHERE) {
	my @thePoses = sort {$relPoses{$b} <=> $relPoses{$a}} keys %relPoses;
	my @absPoses;
	for (@Poses,@thePoses) {
		if ($_ < 0) {
			push @absPoses,-($minLeft - $_);
		} else {
			push @absPoses,($minLeft + $_);
		}
	}
	print "$i->[5]\t$offset\t@{$i}[0..4]; @Poses, @thePoses -> @absPoses\n$seqCIGAR\n";
	print 'B' x ($offset-$firstSC),$i->[9],"\n";
	my $tmpstr = ' ' x length($seqCIGAR);
	for my $p (@Poses) {
		if ($p>=0) {
			substr $tmpstr,$p,1,']';
		} else {
			substr $tmpstr,-$p,1,'[';
		}
	}
	print $tmpstr,"\n";
}
		push @ReadsCIGAR,$seqCIGAR;
	}
	my @thePosesA = sort {$relPoses{$b} <=> $relPoses{$a}} keys %relPoses;
	my $thePos = int($thePosesA[0]);	# $thePos = $thePosesA[0] - 0.5
	my @usingPoses;
	push @usingPoses,$thePos if exists $relPosesFR{$thePos};
	push @usingPoses,-1-$thePos if exists $relPosesFR{-1-$thePos};
	if (@usingPoses == 1) {	# 先找最值处的']'&'['，若没有成对，再在对侧找成对的。只提最近的取一对。
		if ($usingPoses[0] < 0) {
			for my $p ( ($thePos-$main::posAround) .. $thePos ) {
				if (exists $relPosesFR{$p}) {
					push @usingPoses,$p;
				}
			}
			#@usingPoses = @usingPoses[0,-1] if @usingPoses > 2;
			if (@usingPoses > 2) {
				my @t = @usingPoses;
				shift @t;
				@t = sort {$relPosesFR{$b} <=> $relPosesFR{$a} || $b <=> $a } @t;
				@usingPoses = ($usingPoses[0],$t[0]);
			}
		} elsif ($usingPoses[0] > 0) {
			for my $p ( ($thePos+1) .. ($thePos+$main::posAround+1) ) {
				if (exists $relPosesFR{-$p}) {
					push @usingPoses,-$p;
				}
			}
			#@usingPoses = @usingPoses[0,1] if @usingPoses > 2;
			if (@usingPoses > 2) {
				my @t = @usingPoses;
				shift @t;
				@t = sort {$relPosesFR{$b} <=> $relPosesFR{$a} || $a <=> $b } @t;
				@usingPoses = ($usingPoses[0],$t[0]);
			}
		}	# 坐标需要核实下。done.
	}
	my %Bases;
	for my $i (@clipReads) {
		my ($YC) = grep /^YC:Z:/,@$i;
		my $mtype;
		if ($i->[1] & 16) {	# r
			if ($YC eq 'YC:Z:CT') {
				$mtype = 'GA';	# R
			} elsif ($YC eq 'YC:Z:GA') {
				$mtype = 'CT';	# Y
			} else {die;}
		} else {	# f
			if ($YC eq 'YC:Z:CT') {
				$mtype = 'CT';
			} elsif ($YC eq 'YC:Z:GA') {
				$mtype = 'GA';
			} else {die;}
		}
#print "$mtype $i->[9] $YC\n";
		my ($firstSC,$seqCIGAR) = getSeqCIGAR($minLeft,$i);
		my $offset = $i->[3] - $minLeft - $firstSC;
		for my $p (@usingPoses) {
			my ($tlen,$tmp,$vseq,$vqual,$t)=(0);
			if ($p > 0) {	# M]S
				$tmp = substr($seqCIGAR,$p+1) or next;
				# unless ($tmp) {
				# 	ddx $i;
				# }
=pod
11860 sf340_A_Ref_61780248_61780667_61781087_Vir_+_701_820_R_420_90   99      chr1    61780588        60      83M7S   =       61780798        300
11861 sf450_F_Ref_61780248_61780667_61781087_Vir_+_701_820_R_420_90   609     chr1    61780588        0       53S37M  =       61780908        410	Fs90_m458AE.bam，原因不明。MapQ＝0.
11862 sf10_1_Ref_61780248_61780667_61781087_Vir_+_701_820_R_420_90    147     chr1    61780588        60      83M7S   =       61780258        -413
=cut
				$tmp =~ /^(S+)/;
				if (defined $1) {
					$tlen = length $1;
					$t = $p+1 -$offset;
					$vseq = substr $i->[9],$t,$tlen;
					$vqual = substr $i->[10],$t,$tlen;
				}
			} else {	# S[M
				$tmp = substr($seqCIGAR,0,-$p-1) or next;
				$tmp =~ /(S+)$/;
				if (defined $1) {
					$tlen = length $1;
					$t = -$p-1 - $tlen -$offset;
					$vseq = substr $i->[9],$t,$tlen;
					$vqual = substr $i->[10],$t,$tlen;
				}
			}
			if ($vseq and length($vseq) >= $main::minVirLen) {
				push @{$Bases{$p}},[$vseq,$vqual,$mtype];
			}
# print ">>>$tlen, $tmp, $p\n$vseq\n$vqual\n";
# print substr($seqCIGAR,$p,8),"\n" if $p>0;
# print substr($seqCIGAR,-$p-2,8),"\n" if $p<0;
			# }
		}
	}
	ddx \%relPoses,\%relPosesFR,\@usingPoses,\%Bases if $DEBGUHERE;
	my %Results;
	for my $p (@usingPoses) {
		my (@theReads,$absPos);
		if ($p<0) {
			$absPos = $p - $minLeft;
			for (@{$Bases{$p}}) {
				my $x = reverse $_->[0];
				my $y = reverse $_->[1];
				push @theReads,[$x,$y,$_->[2]];
			}
		} else {
			$absPos = $p + $minLeft;
			@theReads = @{$Bases{$p}};
		}
		my $mergstr = mergeStr(\@theReads);
		$mergstr = reverse $mergstr if $p < 0;
		#ddx \@theReads;
		#ddx $Bases{$p};
		my $depth = @theReads;
		#print "$p,$depth\t$mergstr\n";
		$Results{$absPos} = [$depth,$mergstr];
	}
print "@usingPoses\t",'-' x 25,"\n" if $DEBGUHERE;
	# my (%absPoses,%absPosesFR);
	# for (keys %relPoses) {
	# 	$absPoses{$minLeft + $_} = $relPoses{$_};
	# }
	# for (keys %relPosesFR) {
	# 	if ($_ > 0) {
	# 		$absPosesFR{$_ + $minLeft} = $relPosesFR{$_};
	# 	} else {
	# 		$absPosesFR{$_ - $minLeft} = $relPosesFR{$_};
	# 	}
	# }
	return (\%Results);
}
sub mergeStr($) {
	my @Strs = @{$_[0]};
	my ($maxLen,$merged)=(0);
	for (@Strs) {
		my $l = length $_->[0];
		$maxLen = $l if $maxLen < $l;
	}
	for my $p (0 .. ($maxLen-1)) {
		my (%Col,%ColBp) = ();
		for (@Strs) {
			my $str = $_->[0];
			my $qual = $_->[1];
			if ($p < length($str)) {
				my $c = substr $str,$p,1;
				my $q = substr $qual,$p,1;
				if ($_->[2] eq 'CT' and $c eq 'T') {
					$c = 'Y';
				} elsif ($_->[2] eq 'GA' and $c eq 'A') {
					$c = 'R';
				}
				++$ColBp{$c};
				if ($c eq 'Y') {	# CT
					$Col{$_} += $main::Qual2LgP{$q}->[4] for qw(G A);
					$Col{$_} += $main::Qual2LgP{$q}->[3] for qw(C T);
				} elsif ($c eq 'R') {	# GA
					$Col{$_} += $main::Qual2LgP{$q}->[4] for qw(C T);
					$Col{$_} += $main::Qual2LgP{$q}->[3] for qw(G A);
				} else {
					$Col{$_} += $main::Qual2LgP{$q}->[2] for qw(A T C G);
					$Col{$c} += $main::Qual2LgP{$q}->[1];
				}
			}
		}
		my ($res,@Bps);
		if (keys(%ColBp)==1) {
			$res = (keys(%ColBp))[0];
		} else {
			@Bps = sort { $Col{$a} <=> $Col{$b} } keys %Col;	# choose the smaller one.
			if ( $Col{$Bps[1]} - $Col{$Bps[0]} < 0.02 ) {	# q=3时，-10*($lgbp-$lgbq)=0.020624399283,最小。
				my @t = sort { $a cmp $b } @Bps[0,1];
				$res = $REV_IUB{$t[0].$t[1]};
			} else {
				$res = $Bps[0];
			}
		}
		$merged .= $res;
# ddx \%Col,\%ColBp,$res,\@Bps;
	}
	return $merged;
}

sub getDeriv($$$) {
	my ($interest,$bypass,$str) = @_;
	my $len = length $str;
	my $Inserted = 0; # I 与 D 会造成 QSEQ 与 REF 的不同步，必须改成返回两套坐标才行。时间紧迫，先无视。
	my $Deleted = 0; # 现在的 $seqCIGAR 不含‘D’，在遇到D时会左偏。
	my @interestedPoses;	# define as gap after index, which is from 0. Thus cmp between this & next.
	for my $i (0 .. $len-2) {
		my $thischar = substr $str,$i,1;
		my $nextchar = substr $str,$i+1,1;
		#++$Inserted if $thischar eq 'I';
		#++$Deleted if $thischar eq 'D';
		if ($thischar eq $nextchar or $thischar=~/[$bypass]/) {
			next;
		} elsif ($thischar eq $interest) {	# MS -> HV
			push @interestedPoses,$i;
		} elsif ($nextchar eq $interest) {	# SM -> VH
			push @interestedPoses,-($i+1);
		} else {
			next;
		}
	}
	return @interestedPoses;
}

sub mergehash($$) {
	my ($sink,$in) = @_;
	return unless ref($in);
	return unless ref($in) eq "HASH";
	for my $k (keys %{$in}) {
		if (exists $sink->{$k}) {
			$sink->{$k} += $in->{$k};
		} else {
			$sink->{$k} = $in->{$k};
		}
	}
}

1;
