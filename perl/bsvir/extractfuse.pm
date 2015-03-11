use strict;
use warnings;
use File::Basename;
#package main;

# Static Var.
our $DEBUG;
our ($minLen,%Genome,%ChrLen);

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
sub CIAGR2Bin($$$$) {
	my ($seq,$CIGAR,$Value,$acceptdCIAGR)=@_;
	my (@edit_cmd) = $CIGAR =~ m/\d+\w/g;
	my $curr_pos = 0;
	my $total_len =0;
	my $refbin = $seq;
	foreach my $cmd (@edit_cmd) {
		if (my ($M) = $cmd =~ m/(\d+)M/) {
			substr($refbin,$curr_pos,$M,$Value x $M) if $acceptdCIAGR=~/M/;
			print "-M- $refbin,$curr_pos,$total_len,$M,$cmd\n" if $DEBUG > 3;
			$curr_pos += $M;
			$total_len+= $M;
		} elsif (my ($I) = $cmd =~ m/(\d+)I/) {
			if ($acceptdCIAGR=~/I/) {
				substr($refbin,$curr_pos,$I,'');	#delete $I characters
			} else {
				substr($refbin,$curr_pos,$I,$Value x $I);
				$curr_pos += $I;
				$total_len+= $I;
			}
			print "-I- $refbin,$curr_pos,$total_len,$I,$cmd\n" if $DEBUG > 3;
			#$total_len -= $I;
		} elsif (my ($D) = $cmd =~ m/(\d+)D/ and $acceptdCIAGR=~/D/) {
			substr($refbin,$curr_pos,0,$Value x $D);
			print "-D- $refbin,$curr_pos,$total_len,$D,$cmd\n" if $DEBUG > 3;
			$total_len += $D;
			$curr_pos  += $D;
		} elsif (my ($S) = $cmd =~ m/(\d+)S/) {
			print "-S- $refbin,$curr_pos,$total_len,$S,$cmd\n" if $DEBUG > 3;
			# POS │ 1-based leftmost POSition/coordinate of clipped sequence
			# thus no change on $curr_pos and $total_len.
			# It is up to bwa to ensure 'S' exists only at terminals. bsmap does not use 'S'.
			#substr($ref,$curr_pos,$S,'');	#delete $I characters
			substr($refbin,$curr_pos,$S,('0' x $S)) if $acceptdCIAGR=~/s/;
			$curr_pos += $S;
		} elsif ($cmd !~ /\d+[MIDS]/) {die "[x]Unsupported CIAGR:[$cmd]";}
	}
	if ( $DEBUG > 3 ) {
		print "<-$CIGAR,$Value,$acceptdCIAGR\n->$seq\n- $refbin\n";
	}
	return $refbin;
}
sub GenCigar($) {
	my ($seq)=@_;
}

sub Mgfq2HumVir($$) {
	my ($rd,$Type) = @_;	# FQID12_0, FQSeq_1, FQQual, HumFlag_3, HumChr, HumPos, HumCIAGR_6, HumMapQ, VirFlag_8, VirChr, VirPos, VirCIAGR_11, VirMapQ, YDYCHum, YDYCVir_14
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
	my $CIAGRmHum = CIAGR2Bin('0' x length($seqHum),$$rd[3+$idA],'1','MIDS');
	my $CIAGRmVir = CIAGR2Bin('0' x length($seqVir),$$rd[3+$idB],'2','MS');
	$CIAGRmVir = reverse $CIAGRmVir if $strandVir2Hum;
	$CIAGRmVir = CIAGR2Bin($CIAGRmVir,$$rd[3+$idA],'2','IDS');
	my $CIAGRmDiff = $CIAGRmHum ^ $CIAGRmVir;
	$CIAGRmDiff |= '0' x length $CIAGRmDiff;
	GenCigar($CIAGRmDiff);
	if ( $DEBUG > 0 ) {
		print "<<<[$Type]$strandHum,$strandVir2Hum\n ,$$rd[1]\n$strandHum,$seqHum\n$strandVir,$seqVir\n$$rd[6],$$rd[11]\n$CIAGRmHum\n$CIAGRmVir\n$CIAGRmDiff\n";
	}
	return [$strandHum,$strandVir2Hum,$seqHum,$qual];
}

1;
