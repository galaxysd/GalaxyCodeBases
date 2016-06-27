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
sub revcom($) {
	my $str = $_[0];
	$str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
	my $rev = reverse $str;
	$rev =~ tr/[](){}<>/][)(}{></;
	return $rev;
}



sub getsamChrLen($) {
	my $in = $_[0];
	my $inpath = dirname($in);
	my ($Ref,$fq1,$fq2,$readlen1,$readlen2);
	open( IN,"-|","samtools view -H $in") or die "Error(m1) opening $in: $!\n";
	while (<IN>) {
		chomp;
		my ($t,$id,$ln)=split /\t/;
		if ($t eq '@SQ') {
			$id = (split /\:/,$id)[1];
			$ln = (split /\:/,$ln)[1];
			$ChrLen{$id} = $ln;
		} elsif ($t eq '@PG') {
# @PG	ID:BSMAP	VN:2.87	CL:"bsmap -u -d HomoGRCh38.fa -a F12HPCCCSZ0010_Upload/s00_C.bs_1.fq.gz -b F12HPCCCSZ0010_Upload/s00_C.bs_2.fq.gz -z 64 -p 12 -v 10 -q 2 -o s00_C.bs.bam"
# @PG	ID:bwa-meth	PN:bwa-meth	VN:0.10	CL:"./bwameth.py -t 24 -p s00_C.bshum --read-group s00_C --reference HomoGRCh38.fa s00_C.bs_1.fq.gz s00_C.bs_2.fq.gz"
			if ($id eq 'ID:BSMAP') {
				$Ref = $1 if /\-d\s+([.\w]+)\b/;
				$fq1 = $1 if /\-a\s+([^\s"]+)(\s|"?$)/;
				$fq2 = $1 if /\-b\s+([^\s"]+)(\s|"?$)/;
			} elsif ($id eq 'ID:bwa-meth') {
				$Ref = $1 if /\--reference\s+([.\w]+)\b/;
				if (/CL:.+\s([^-\s"]+)\s+([^-\s"]+)(\s|"?$)/) {
					$fq1 = $1;
					$fq2 = $2;
				}
			}

		}
	}
	if ($inpath ne '.') {
		($Ref,$fq1,$fq2) = map { "$inpath/$_" } ($Ref,$fq1,$fq2);
	}
#warn "$Ref,$fq1,$fq2,$inpath";
	my $FQ1 = openfile($fq1);
	<$FQ1>;chomp($_=<$FQ1>);
	$readlen1 = length $_;
	close $FQ1;
	my $FQ2 = openfile($fq2);
	<$FQ2>;chomp($_=<$FQ2>);
	$readlen2 = length $_;
	close $FQ2;
	if ( ($DEBUG - int($DEBUG)) < 0.85 or ($Ref =~ /HBV\.AJ507799\.2\.fa$/) ) {
		my $GENOME = openfile($Ref);
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
			if ( $DEBUG > 0 ) {
				die if $thelength != $ChrLen{$seqname};
			}
		}
		close $GENOME;
	}
	return [$Ref,$fq1,$fq2,$readlen1,$readlen2];
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
		}
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



1;
