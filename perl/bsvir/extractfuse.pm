use strict;
use warnings;
use File::Basename;
#package main;

# Static Var.
our $DEBUG;
our ($minLen,%Genome,%ChrLen);
our ($low_complexity_cutoff_ratio,$N_num_cutoff);

# Functions
sub openfile($) {
	my ($filename)=@_;
	my $infile;
	if ($filename=~/.xz$/) {
		open( $infile,"-|","xz -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.gz$/) {
		open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.bz2$/) {
		open( $infile,"-|","bzip2 -dc $filename") or die "Error opening $filename: $!\n";
	} else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
	return $infile;
}

sub getRefChrLen($) {
	my $in = $_[0];
	if ( ($DEBUG - int($DEBUG)) < 0.85 or ($in =~ /HBV\.AJ507799\.2\.fa$/) ) {
		my $GENOME = openfile($in);
		while (<$GENOME>) {
			s/^>//;
			/^(\S+)/ or next;
			my $seqname = $1;
			print STDERR " >$seqname ...";
			$/=">";
			my $genome=<$GENOME>;
			chomp $genome;
			$genome=~s/\s//g;
			$/="\n";
			$Genome{$seqname}=$genome;
			my $thelength = length $Genome{$seqname};
			print STDERR "\b\b\b", $thelength, ".\n";
			$genome='';
			$ChrLen{$seqname} = $thelength;
		}
		close $GENOME;
	}
}

sub revcom($) {
	my $str = $_[0];
	$str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
	my $rev = reverse $str;
	$rev =~ tr/[](){}<>/][)(}{></;
	return $rev;
}

sub parseCIGAR($$) {
	# http://davetang.org/muse/2011/01/28/perl-and-sam/	<= insert '-' to both
	# http://www.perlmonks.org/?node_id=858230	<= insert X when 'D' in reads. <- chosen
	my ($Read,$CIGAR)=@_;
	my $ref = $Read;
	my (@edit_cmd) = $CIGAR =~ m/\d+\w/g;
	my $curr_pos = 0;
	my $total_len =0;
	foreach my $cmd (@edit_cmd) {
		#warn "$cmd";
		if (my ($M) = $cmd =~ m/(\d+)M/) {
			$curr_pos += $M;
			$total_len+= $M;
		} elsif (my ($I) = $cmd =~ m/(\d+)I/) {
			substr($ref,$curr_pos,$I,'');	#delete $I characters
			#$total_len -= $I;
		} elsif (my ($D) = $cmd =~ m/(\d+)D/) {
			substr($ref,$curr_pos,0,"!" x $D);  #insert $D as '!', which is Q=0 for QUAL
			$total_len += $D;
			$curr_pos  += $D;
		} elsif (my ($S) = $cmd =~ m/(\d+)S/) {
			# POS │ 1-based leftmost POSition/coordinate of clipped sequence
			# thus no change on $curr_pos and $total_len.
			# It is up to bwa to ensure 'S' exists only at terminals. bsmap does not use 'S'.
			substr($ref,$curr_pos,$S,'');	#delete $I characters
			#substr($ref,$curr_pos,$S,(' ' x $S));
			#$curr_pos += $S;
		} else {die;}
	}
	#$ref =~ s/ //g;
	#$ref = substr($ref,0,$total_len); #truncate ?????
	if ( $DEBUG > 0 ) {
		if (length $ref != $total_len) {
			print "INPUT = [$Read], CIGAR=[$CIGAR]\n";
			print "->REF = [$ref]\n";
			print length($ref),"->$total_len\n\n";
			die;
		}
	}
	return $ref;
}
sub CIGAR2Bin($$$$) {
	my ($seq,$CIGAR,$Value,$acceptdCIGAR)=@_;
	my (@edit_cmd) = $CIGAR =~ m/\d+\w/g;
	my $curr_pos = 0;
	my $total_len =0;
	my $refbin = $seq;
	foreach my $cmd (@edit_cmd) {
		if (my ($M) = $cmd =~ m/(\d+)M/) {
			substr($refbin,$curr_pos,$M,$Value x $M) if $acceptdCIGAR=~/M/;
			print "-M- $refbin,$curr_pos,$total_len,$M,$cmd\n" if $DEBUG > 3;
			$curr_pos += $M;
			$total_len+= $M;
		} elsif (my ($I) = $cmd =~ m/(\d+)I/) {
			if ($acceptdCIGAR=~/I/) {
				substr($refbin,$curr_pos,$I,'');	#delete $I characters
			} else {
				substr($refbin,$curr_pos,$I,$Value x $I);
				$curr_pos += $I;
				$total_len+= $I;
			}
			print "-I- $refbin,$curr_pos,$total_len,$I,$cmd\n" if $DEBUG > 3;
			#$total_len -= $I;
		} elsif (my ($D) = $cmd =~ m/(\d+)D/ and $acceptdCIGAR=~/D/) {
			my $thisValue = $Value;
			if ($Value == -1) {
				$thisValue = substr($refbin,$curr_pos,1);
			}
			substr($refbin,$curr_pos,0,$thisValue x $D);
			print "-D- $refbin,$curr_pos,$total_len,$D,$cmd,[$thisValue]\n" if $DEBUG > 3;
			$total_len += $D;
			$curr_pos  += $D;
		} elsif (my ($S) = $cmd =~ m/(\d+)S/) {
			print "-S- $refbin,$curr_pos,$total_len,$S,$cmd\n" if $DEBUG > 3;
			# POS │ 1-based leftmost POSition/coordinate of clipped sequence
			# thus no change on $curr_pos and $total_len.
			# It is up to bwa to ensure 'S' exists only at terminals. bsmap does not use 'S'.
			#substr($ref,$curr_pos,$S,'');	#delete $I characters
			substr($refbin,$curr_pos,$S,('0' x $S)) if $acceptdCIGAR=~/s/;
			$curr_pos += $S;
		} elsif ($cmd !~ /\d+[MIDS]/) {die "[x]Unsupported CIGAR:[$cmd]";}
	}
	if ( $DEBUG > 3 ) {
		print "<-$CIGAR,$Value,$acceptdCIGAR\n->$seq\n- $refbin\n";
	}
	return $refbin;
}
sub GenRLEa($) {
	my ($seq)=@_;
	my (@retChar,@retLen,@retLastP);
	my $lastPos = 0;
	while ($seq =~ m/((\d)\2*)/g) {
		#$ret .= $bin2ab{$2} . length($1);
		push @retChar,"$2";
		push @retLen,length($1);
		push @retLastP,$lastPos;
		$lastPos += length($1);
	}
	return [\@retChar,\@retLen,\@retLastP];
}
sub RLEa2str($) {
	my $dat = $_[0];
	my ($retCharR,$retLenR) = @$dat;
	my %bin2ab = (0=>'N',1=>'A',2=>'B',3=>'C');
	my $ret;
	for my $i (0 .. $#$retCharR) {
		if (exists $bin2ab{$$retCharR[$i]}) {
			$ret .= $bin2ab{$$retCharR[$i]} . $$retLenR[$i];
		} else {die;}
	}
	return $ret;
}
sub low_complexity_filter($) {  # https://github.com/chienchi/FaQCs/blob/master/FaQCs.pl
	my ($seq) = @_;
	my $seq_len = length ($seq);
	my @low_complex_array=("A","T","C","G","AT","AC","AG","TA","TC","TG","CA","CT","CG","GA","GT","GC","N");
	my $low_complexity=0;
	for my $item (@low_complex_array) {
		my $item_len = length ($item);
		my $num_low_complex = $seq =~ s/$item/$item/g;
		my $value = ($num_low_complex*$item_len/$seq_len);
		if ( $value >= $low_complexity) {
			$low_complexity=$value;
			last if $value == 1;
		}
	}
	return ($low_complexity);
}
sub analyseRLEa($$) {
	my ($dat,$seq)=@_;
	my ($retCharR,$retLenR,$retLastPR) = @$dat;
	my @ret;
	if ($#$retCharR==2 and $$retCharR[1]==3) {
		if (($$retCharR[0]==1 and $$retCharR[2]==2) or ($$retCharR[0]==2 and $$retCharR[2]==1)) {
			my $subseq = substr $seq,$$retLastPR[1],$$retLenR[1];
			my $low_complexity = low_complexity_filter($subseq);
			my $drop=0;
			$drop = 1 if $low_complexity > $low_complexity_cutoff_ratio;
			if ( $DEBUG > 2 ) {
				print "-!- $seq,",RLEa2str($dat),"\n-!- ",'.'x$$retLenR[0],"$subseq,$drop,$low_complexity\n";
			}
			unless ($drop) {
				@ret = ($$retLastPR[1],$$retLastPR[2]-1);
			}
		}
	}
	return \@ret;
}
sub Mgfq2HumVir($$) {
	my ($rd,$Type) = @_;	# FQID12_0, FQSeq_1, FQQual, HumFlag_3, HumChr, HumPos, HumCIGAR_6, HumMapQ, VirFlag_8, VirChr, VirPos, VirCIGAR_11, VirMapQ, YDYCHum, YDYCVir_14
	my ($strandHum,$strandVir,$idA,$idB)=(0,0);
	if ($Type eq 'Hum') {
		($idA,$idB) = (3,8);	# Hum is 'Hum'
	} elsif ($Type eq 'Vir') {
		($idA,$idB) = (8,3);	# Hum is 'Vir'
	} else {die;}
	$strandHum = 1 if $$rd[$idA] & 0x10;
	$strandVir = 1 if $$rd[$idB] & 0x10;
	my $strandVir2Hum = ($strandVir ^ $strandHum);
	my ($seqHum,$seqVir,$qual)=($$rd[1],$$rd[1],$$rd[2]);
	if ($strandVir or $strandHum) {
		my $rcseq = revcom($$rd[1]);
		if ($strandHum) {
			$seqHum = $rcseq;
			$qual = reverse $$rd[2];
		}
		$seqVir = $rcseq if $strandVir;
	}
	my $CIGARmHum = CIGAR2Bin('0' x length($seqHum),$$rd[3+$idA],'1','MIDS');
	my $CIGARmVir = CIGAR2Bin('0' x length($seqVir),$$rd[3+$idB],'2','MS');
	$CIGARmVir = reverse $CIGARmVir if $strandVir2Hum;
	$CIGARmVir = CIGAR2Bin($CIGARmVir,$$rd[3+$idA],-1,'IDS');
	my $CIGARmDiff = $CIGARmHum ^ $CIGARmVir;
	$CIGARmDiff |= '0' x length $CIGARmDiff;
	my $retRLEa = GenRLEa($CIGARmDiff);
	my $retana = analyseRLEa($retRLEa,$seqHum);
	if ( $DEBUG > 2 or ($DEBUG > 1 and $#$retana>-1) ) {
		print "<<<[$Type]$strandHum,$strandVir2Hum\n ,$$rd[1]\n$strandHum,$seqHum\n$strandVir,$seqVir\n$$rd[6],$$rd[11]\n$CIGARmHum\n$CIGARmVir\n$CIGARmDiff,",
		RLEa2str($retRLEa),"\n";
		ddx $retRLEa;
		ddx $retana;
	}
	return [$retana,$strandHum,$strandVir2Hum,$retRLEa,$seqHum,$qual];
}

1;
