#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

my $minLen = 5;
my $DEBUG=2.9;

die "Usage: $0 <ins_size> <hum_sorted_bam> <vir_sorted_bam> [out_prefix]\n" if @ARGV < 3;
my ($isize,$inhum,$invir,$out)=@ARGV;

unless (defined $out) {
	$out='fuse_'.$inhum;
	$out =~ s/\.[sb]am(\.gz)?//g;
}
warn "[!]IN: h=[$inhum],v=[$invir],i=[$isize] to [$out].*\n";

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

my %Genome;	# Let's share it.
sub getsamChrLen($) {
	my $in = $_[0];
	my (%ChrLen,$Ref,$fq1,$fq2,$readlen1,$readlen2);
	open( IN,"-|","samtools view -H $in") or die "Error opening $in: $!\n";
	while (<IN>) {
		chomp;
		my ($t,$id,$ln)=split /\t/;
		if ($t eq '@SQ') {
			$id = (split /\:/,$id)[1];
			$ln = (split /\:/,$ln)[1];
			$ChrLen{$id} = $ln;
		} elsif ($t eq '@PG') {
# @PG	ID:BSMAP	VN:2.87	CL:"bsmap -u -d HomoGRCh38.fa -a F12HPCCCSZ0010_Upload/s00_C.bs_1.fq.gz -b F12HPCCCSZ0010_Upload/s00_C.bs_2.fq.gz -z 64 -p 12 -v 10 -q 2 -o s00_C.bs.bam"
			$Ref = $1 if /\-d\s+([.\w]+)\b/;
			$fq1 = $1 if /\-a\s+([^\s"]+)(\s|"?$)/;
			$fq2 = $1 if /\-b\s+([^\s"]+)(\s|"?$)/;
		}
	}
	my $FQ1 = openfile($fq1);
	<$FQ1>;chomp($_=<$FQ1>);
	$readlen1 = length $_;
	close $FQ1;
	my $FQ2 = openfile($fq2);
	<$FQ2>;chomp($_=<$FQ2>);
	$readlen2 = length $_;
	close $FQ2;
	unless (($DEBUG - int($DEBUG)) > 0.8) {
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
	return [[$Ref,$fq1,$fq2,$readlen1,$readlen2],%ChrLen];
}

my ($ret,%ChrLenHum,$finHum,%ChrLenHBV,$finHBV);

$ret = getsamChrLen($inhum);
$finHum = shift @$ret; %ChrLenHum = @$ret;

$ret = getsamChrLen($invir);
$finHBV = shift @$ret; %ChrLenHBV = @$ret;

#ddx \%ChrLenHum;
ddx \%ChrLenHBV;
ddx $finHum;
ddx $finHBV;

warn "[!]Ref: h=[$$finHum[0]],v=[$$finHBV[0]]\n";

sub parseCIGAR($$) {
	# http://davetang.org/muse/2011/01/28/perl-and-sam/	<= insert '-' to both
	# http://www.perlmonks.org/?node_id=858230	<= insert X when 'D' in reads. <- chosen
	my ($Read,$CIGAR)=@_;
	#my ($input, $ref) = ($Read, $Read);
	my $ref = $Read;
	my (@edit_cmd) = $CIGAR =~ m/\d+\w/g;
	my $curr_pos = 0;
	my $total_len =0;
	foreach my $cmd (@edit_cmd) {
		if (my ($M) = $cmd =~ m/(\d+)M/) {
			$curr_pos += $M;
			$total_len+= $M;
		} elsif (my ($I) = $CIGAR =~ m/(\d+)I/) {
			substr($ref,$curr_pos,$I,'');	#delete $I characters
			$total_len -= $I;
		} elsif (my ($D) = $CIGAR =~ m/(\d+)D/) {
			substr($ref,$curr_pos,0,"!" x $D);  #insert $D as '!', which is Q=0 for QUAL
			$total_len += $D;
			$curr_pos  += $D;
		} elsif (my ($S) = $CIGAR =~ m/(\d+)S/) {
			# POS │ 1-based leftmost POSition/coordinate of clipped sequence
			# thus no change on $curr_pos and $total_len.
			# It is up to bwa to ensure 'S' exists only at terminals. bsmap does not use 'S'.
			substr($ref,$curr_pos,$S,'');	#delete $I characters
			print STDERR '.1.';
		}
	}
	#$ref = substr($ref,0,$total_len); #truncate ?????
	if ( $DEBUG > 0 ) {
		die if length $ref != $total_len;
	}
	#print "INPUT = $input   CIGAR = $CIGAR\n";
	#print "->REF = $ref\n\n";
	return $ref;
}

open INVIR,"-|","samtools view $invir" or die "Error opening $invir: $!\n";
my (%PileUpVIR,%DepthVIR);
while (<INVIR>) {
	chomp;
	my ($id, $flag, $ref, $pos, $mapq, $CIAGR, $mref, $mpos, $isize, $seq, $qual, @OPT) = split /\t/;
	my $mappedSeq = parseCIGAR($seq,$CIAGR);
	my $mappedQual = parseCIGAR($qual,$CIAGR);
	my $mappedLen = length $mappedSeq;
	if ( $DEBUG > 1 ) {
		push @{$PileUpVIR{$ref}},[$pos,$mappedLen,$id,$mappedSeq,$mappedQual,$seq,$CIAGR];
	} else {
		push @{$PileUpVIR{$ref}},[$pos,$mappedLen,$id,$mappedSeq,$mappedQual];
	}
}
for my $chr (keys %PileUpVIR) {
	for my $Datp ( @{$PileUpVIR{$chr}} ) {
		for my $pos ( ($Datp->[0]) .. ($Datp->[0] + $Datp->[1]) ) {
			++$DepthVIR{$chr}->[$pos];
		}
	}
}

ddx \%DepthVIR; 

__END__
./getfuse.pl 200 grep_s00_C.bs.virsam.sort.bam grep_s00_C.bs.sort.bam
记得把流程中的文件名改下，突出ref是哪个物种。


CMD Log:
bsmap -u -d HomoGRCh38.fa -a F12HPCCCSZ0010_Upload/s00_C.bs_1.fq.gz -b F12HPCCCSZ0010_Upload/s00_C.bs_2.fq.gz -z 64 -p 12 -v 10 -q 2 -o s00_C.bs.bam >s00_C.bs.log 2>s00_C.bs.err

./getUnPaired.pl s00_C.bs.bam	# 出UnPaired的fq和sam，及grep_s00_C.insertsize。

bsmap -z 64 -p 12 -v 10 -q 2 -d HBV.AJ507799.2.fa -a grep_s00_C.bs.1.fq.gz -b grep_s00_C.bs.2.fq.gz -o grep_s00_C.bs.bam 2>grep_s00_C.bs.log &
samtools sort -l 9 grep_s00_C.bs.bam grep_s00_C.bs.sort	# UnPaired与病毒比

samtools view -bS grep_s00_C.bs.virsam.gz > grep_s00_C.bs.virsam.bam && \
samtools sort -l 9 grep_s00_C.bs.virsam.bam grep_s00_C.bs.virsam.sort &	# UnPaired与人比

./getfuse.pl 200 grep_s00_C.bs.virsam.sort.bam grep_s00_C.bs.sort.bam
