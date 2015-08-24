package main;
use strict;
use warnings;

sub revcom($) {
	my $str = $_[0];
	$str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
	my $rev = reverse $str;
	$rev =~ tr/[](){}<>/][)(}{></;
	return $rev;
}

sub openfile($) {
	my ($filename)=@_;
	my $infile;
	if ($filename=~/.[bs]am(.gz)?$/) {
		open( $infile,"-|","samtools view -F 256 $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.gz$/) {
		open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.bz2$/) {
		open( $infile,"-|","bzip2 -dc $filename") or die "Error opening $filename: $!\n";
	} else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
	return $infile;
}

sub getRefChr1stID($) {
	my $GENOME = $_[0];
	while (<$GENOME>) {
		s/^>//;
		/^(\S+)/ or next;
		my $seqname = $1;
		return $seqname;
	}
}

sub getRefChr1st($) {
	my $GENOME = $_[0];
	while (<$GENOME>) {
		s/^>//;
		/^(\S+)/ or next;
		my $seqname = $1;
		print STDERR ">$seqname ...";
		$/=">";
		my $genome=<$GENOME>;
		chomp $genome;
		$genome=~s/\s//g;
		$/="\n";
		my $thelength = length $genome;
		print STDERR "\b\b\b", $thelength, ".\n";
		return $genome;
	}
}

sub getticks($$$$$) {
	my ($RefBorder,$Refstr,$RefLen,$PEinsertLen,$RefNratioMax) = @_;
	my @theticks;
	while (@theticks < 100) {
		my $pos0 = int(rand($RefLen-(2*$RefBorder)))+$RefBorder;
		my $str0 = substr $Refstr,($pos0-$PEinsertLen),2*$PEinsertLen;
		my $seq = $str0;
		my $N = $seq=~tr/Nn//;
		next if $N > 2*$PEinsertLen*$RefNratioMax;
		push @theticks,$pos0
	}
	@theticks = sort {$a<=>$b} @theticks;
	#warn "$RefBorder,$RefLen";
	return \@theticks;
}

sub getype($$) {
	my ($R1,$R2)=@_;
	my $type='I'; # !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ
	if ($R1 =~ /^[A-Z]+$/) {
		if ($R2 =~ /^[A-Z]+[a-z]+$/) {
			$type = '1';
		} elsif ($R2 =~ /^[A-Z]+[a-z]+[A-Z]+$/) {
			$type = '2';
		} elsif ($R2 =~ /^[a-z]+$/) {
			$type = '3';
		} elsif ($R2 =~ /^[a-z]+[A-Z]+$/) {
			$type = '4';
		} elsif ($R2 =~ /^[A-Z]+$/) {
			$type = '5';
		}
	} elsif ($R1 =~ /^[A-Z]+[a-z]+$/) {
		if ($R2 =~ /^[A-Z]+[a-z]+$/) {
			$type = '6';
		} elsif ($R2 =~ /^[A-Z]+[a-z]+[A-Z]+$/) {
			$type = '7';
		} elsif ($R2 =~ /^[a-z]+$/) {
			$type = '8';
		} elsif ($R2 =~ /^[a-z]+[A-Z]+$/) {
			$type = '9';
		} elsif ($R2 =~ /^[A-Z]+$/) {
			$type = 'A';
		}
	} elsif ($R1 =~ /^[A-Z]+[a-z]+[A-Z]+$/) {
		if ($R2 =~ /^[A-Z]+[a-z]+[A-Z]+$/) {
			$type = 'B';
		} elsif ($R2 =~ /^[a-z]+[A-Z]+$/) {
			$type = 'C';
		} elsif ($R2 =~ /^[A-Z]+$/) {
			$type = 'D';
		}
	} elsif ($R1 =~ /^[a-z]+$/) {
		if ($R2 =~ /^[a-z]+[A-Z]+$/) {
			$type = 'E';
		} elsif ($R2 =~ /^[A-Z]+$/) {
			$type = 'F';
		}
	} elsif ($R1 =~ /^[a-z]+[A-Z]+$/) {
		if ($R2 =~ /^[a-z]+[A-Z]+$/) {
			$type = 'H';
		} elsif ($R2 =~ /^[A-Z]+$/) {
			$type = 'G';
		}
	}
	return $type;
}

sub dosim($$$) {
	my ($Refstr,$Virstr,$Paras)=@_;
	my $PEinsertLen = $Paras->{PEinsertLen};
	my $SeqReadLen = $Paras->{SeqReadLen};
	my $RefBorder = $PEinsertLen + 1000;
	open O,'>',$Paras->{OutPrefix}.'.Ref.fa';
	open R1,'>',$Paras->{OutPrefix}.'.1.fq';
	open R2,'>',$Paras->{OutPrefix}.'.2.fq';
	my @Refticks = @{getticks($RefBorder,$Refstr,$Paras->{RefLen},$PEinsertLen,$Paras->{RefNratioMax})};
	my @Virticks = @{getticks($Paras->{VirFrag},$Virstr,$Paras->{VirLen},int(0.9+ 0.5*$Paras->{VirFrag}),$Paras->{RefNratioMax})};
	ddx $Paras;
	for my $pRef (@Refticks) {
		my $seqR1 = substr $Refstr,($pRef-$PEinsertLen),$PEinsertLen;
		my $seqR2 = substr $Refstr,$pRef,$PEinsertLen;
		my $pVir = shift @Virticks;
		my $startV = $pVir-int(0.5*$Paras->{VirFrag});
		my $seqV = substr $Virstr,$startV,$Paras->{VirFrag};
		my $isReverse = int(rand(2));
		my $strand = '+';
		if ($isReverse) {
			$seqV = revcom($seqV);
			$strand = '-';
		}
		my $newSeq = join('',uc $seqR1,lc $seqV,uc $seqR2);
		my $tID = join('_','Ref',$pRef-$PEinsertLen+1,$pRef,$pRef+$PEinsertLen,'Vir',$strand,$startV+1,$startV+$Paras->{VirFrag},'R',$PEinsertLen,$SeqReadLen);
		print O '>',$tID,"\n$newSeq\n\n";
		my $maxP = length($newSeq) - $PEinsertLen;
		#for my $p ($Paras->{LeftStart} .. $Paras->{LeftEnd}) {
		for my $p (0 .. $maxP) {
			#last if $p > $maxP;
			my $PE = substr $newSeq,$p,$PEinsertLen;
			my $R1 = substr $PE,0,$SeqReadLen;
			my $R2 = substr $PE,$PEinsertLen-$SeqReadLen,$SeqReadLen;
			my $revR2 = revcom($R2);
			my $type = getype($R1,$R2);
			$type = '0' if $p == 0 or $p == $maxP;
			my ($Part1,$Part2);
			my $Qual = $type x $SeqReadLen;
			$Part1 = join '-',getInsertParts($PEinsertLen,$SeqReadLen,$Paras->{VirFrag},'f',$p,1);
			$Part2 = join '-',getInsertParts($PEinsertLen,$SeqReadLen,$Paras->{VirFrag},'f',$p,2);
			print R1 "\@sf${p}_${tID}/1 ${Part1} $type\n$R1\n+\n$Qual\n";
			print R2 "\@sf${p}_${tID}/2 ${Part2} $type\n$revR2\n+\n$Qual\n";
			my $revR1 = revcom($R1);
			$type =~ tr/123456789ABCDEFGH/GDFA5HCE94B728316/;	# 反向后的对应关系
			$Qual = $type x $SeqReadLen;
			$Part2 = join '-',getInsertParts($PEinsertLen,$SeqReadLen,$Paras->{VirFrag},'r',$p,2);
			$Part1 = join '-',getInsertParts($PEinsertLen,$SeqReadLen,$Paras->{VirFrag},'r',$p,1);
			print R2 "\@sr${p}_${tID}/2 ${Part2} $type\n$revR1\n+\n$Qual\n";
			print R1 "\@sr${p}_${tID}/1 ${Part1} $type\n$R2\n+\n$Qual\n";
		}
	}
	close O;
	close R1; close R2;
}

sub cigar2rpos($$) {
	my ($cigar,$lpos) = @_;
	my $reflen = 0;
	my @cigar = $cigar =~ /(\d+)(\w)/g;
	while (@cigar) {
		my ($len,$op) = splice(@cigar,0,2);
		$reflen += $len if $op =~ /^[MD]$/;
	}
	return $lpos + $reflen -1;
}

sub cigar2SMS($) {
	my ($cigar) = @_;
	my @cigars = $cigar =~ /(\d+\w)/g;
	my ($Left,$Middle,$Right)=(0,0,0);
	if ($cigars[0] =~ /(\d+)S/) {
		shift @cigars;
		$Left = $1;
		
	}
	if ($cigars[-1] =~ /(\d+)S/) {
		pop @cigars;
		$Right = $1;
	}
	$Middle = join('',@cigars);
	$Middle = cigar2rpos($Middle,1);
	return [$Left,$Middle,$Right];
}



sub getInsertPos($$$$$) {	# 返回Read在模拟拼合片段上的左右端点。
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
sub InsertPos2PartLVR($$$) {	# 对getInsertPos返回的，模拟拼合片段上的点，返回其所属区段:[LVR],及到该区段末尾的长度（与读长无关）。
	my ($Pos,$InsertSize,$VirFrag)=@_;
	my ($type,$lastingLen)=('',0);

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
=pod
	if ($r1fr eq 'f') {
		# those above
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
=cut
	return ($type,$lastingLen);
}
sub Parts2List($$$$$$$$) {	# 根据左右两个InsertPos2PartLVR返回值，计算理论覆盖模式。 [mmL,nnV,ddR] (mm+nn+dd = ReadLength)
	my ($Pos,$type1,$lastingLen1,$type2,$lastingLen2,$ReadLen,$InsertSize,$VirFrag) = @_;
	my %Type2Int = (
		L => 1,
		V => 2,
		R => 3,
	);
	my %Int2Type = reverse %Type2Int;
	#if ($type1 eq $type2) {	# 'L'的Read1记录这样到LR分界点最方便。需要`$nextLL`则手动替换。
	#	return ( "${ReadLen}${type1}" );
	#} els
	if ( $Type2Int{$type1} <= $Type2Int{$type2} ) {
		my @ret = ( "${lastingLen1}${type1}" );
		my $nextInt = 1+ $Type2Int{$type1};
		my $nextPos = $Pos;
		my $nextLL = $lastingLen1;
		my $nextType;
		my $RemainLen = $ReadLen;
		while ($nextInt <= $Type2Int{$type2}) {
			$nextPos += $nextLL;
			$RemainLen -= $nextLL;
			($nextType,$nextLL) = InsertPos2PartLVR($nextPos,$InsertSize,$VirFrag);
			$nextLL=$RemainLen if $nextInt == $Type2Int{$type2};
			#ddx [$nextInt,$nextType,$RemainLen,$nextLL,1,$r1fr,$Pos,$type1,$lastingLen1,$type2,$lastingLen2,$ReadLen,$InsertSize,$VirFrag];
			push @ret,"${nextLL}$nextType";
			++$nextInt;
		}
		return @ret;
	} else { die 'E'; }
}
sub InsertPos2InsertParts($$$$$) {
	my ($InsertSize,$ReadLen,$VirFrag,$FiveT,$ThreeT)=@_;
	my ($type5,$lastingLen5) = InsertPos2PartLVR($FiveT,$InsertSize,$VirFrag);
	my ($type3,$lastingLen3) = InsertPos2PartLVR($ThreeT,$InsertSize,$VirFrag);
	my @Parts = Parts2List($FiveT,$type5,$lastingLen5,$type3,$lastingLen3,$ReadLen,$InsertSize,$VirFrag);
	return @Parts;
}
sub getInsertParts($$$$$$) {
	my ($InsertSize,$ReadLen,$VirFrag,$r1fr,$innerPos,$r12) = @_;
	my ($FiveT,$ThreeT) = getInsertPos($r1fr,$innerPos,$InsertSize,$ReadLen,$r12);
	my @Parts = InsertPos2InsertParts($InsertSize,$ReadLen,$VirFrag, $FiveT,$ThreeT);
	return @Parts;
}


1;
