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

sub dosim($$$) {
	my ($Refstr,$Virstr,$Paras)=@_;
	my $PEinsertLen = $Paras->{PEinsertLen};
	my $SeqReadLen = $Paras->{SeqReadLen};
	open O,'>',$Paras->{OutPrefix}.'.Ref.fa';
	open R1,'>',$Paras->{OutPrefix}.'.1.fq';
	open R2,'>',$Paras->{OutPrefix}.'.2.fq';
	my @Refticks = @{getticks($Paras->{RefBorder},$Refstr,$Paras->{RefLen},$PEinsertLen,$Paras->{RefNratioMax})};
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
		my $tID = join('_','Ref',$pRef-$PEinsertLen,$pRef,$pRef+$PEinsertLen,'Vir',$strand,$startV,$startV+$Paras->{VirFrag});
		print O '>',$tID,"\n$newSeq\n\n";
		my $maxP = length($newSeq) - $PEinsertLen;
		for my $p ($Paras->{LeftStart} .. $Paras->{LeftEnd}) {
			last if $p > $maxP;
			my $PE = substr $newSeq,$p,$PEinsertLen;
			my $R1 = substr $PE,0,$SeqReadLen;
			my $R2 = substr $PE,$PEinsertLen-$SeqReadLen,$SeqReadLen;
			#my $revR1 = revcom($R1);
			my $revR2 = revcom($R2);
			my $Qual = 'e' x $SeqReadLen;
			print R1 "\@sf${p}_${tID}/1\n$R1\n+\n$Qual\n";
			print R2 "\@sf${p}_${tID}/2\n$revR2\n+\n$Qual\n";
			#print R2 "\@sr${p}_${tID}/2\n$revR1\n+\n$Qual\n";
			#print R1 "\@sr${p}_${tID}/1\n$R2\n+\n$Qual\n";
		}
	}
	close O;
	close R1; close R2;
}

1;
