#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;

my $minSoftClipLen = 15;

die "Usage: $0 <in_bam> [out_prefix] [fq1] [fq2]\n" if @ARGV < 1;
my ($in,$out,$fq1,$fq2)=@ARGV;
unless (defined $out) {
	$out='grep_'.$in;
	$out =~ s/\.[sb]am(\.gz)?//g;
}

sub openpipe($$) {
	my ($cmd,$filename)=@_;
	my $infile;	# undef
	if ( length($filename) ) {
		open( $infile,"|-","$cmd >$filename") or die "Error opening [$cmd],[$filename]: $!\n";
	}
	return $infile;
}
my $OUT0 = openpipe('gzip -9c',"$out.vircandi.sam.gz");
my $OUT4 = openpipe('gzip -9c',"$out.shortsoftclip.sam.gz");
open OUT3,'>',"$out.insertsize" or die "Error opening [$out.insertsize]: $!\n";

my $inType='BSMAP';

# unless (defined $fq2) {
	open( IN,"-|","samtools view -H $in") or die "Error opening $in: $!\n";
	my $inpath = dirname($in);
	while (my $line = <IN>) {
		print $OUT0 $line;
		print $OUT4 $line;
		if ( (! defined $fq2) and $line =~ /^\@PG/) {
			$line =~ /CL:"([^"]+)"/ or die "[x]SAM/BAM file header error as [$line]\n";
			my ($t,$id)=split /\t/,$line;
# @PG	ID:BSMAP	VN:2.87	CL:"bsmap -u -d HomoGRCh38.fa -a F12HPCCCSZ0010_Upload/s00_C.bs_1.fq.gz -b F12HPCCCSZ0010_Upload/s00_C.bs_2.fq.gz -z 64 -p 12 -v 10 -q 2 -o s00_C.bs.bam"
# @PG	ID:bwa-meth	PN:bwa-meth	VN:0.10	CL:"./bwameth.py -t 24 -p s00_C.bshum --read-group s00_C --reference HomoGRCh38.fa s00_C.bs_1.fq.gz s00_C.bs_2.fq.gz"
			if ($id eq 'ID:BSMAP') {
				$fq1 = $1 if $line =~ /\-a\s+([^\s"]+)(\s|"?$)/;
				$fq2 = $1 if $line =~ /\-b\s+([^\s"]+)(\s|"?$)/;
			} elsif ($id eq 'ID:bwa-meth') {
				$inType='bwa-meth';
				if ($line =~ /CL:.+\s([^-\s"]+)\s+([^-\s"]+)(\s|"?$)/) {
					$fq1 = $1;
					$fq2 = $2;
				}
			}
			if ($inpath ne '.') {
				($fq1,$fq2) = map { "$inpath/$_" } ($fq1,$fq2);
			}
		}
	}
	close IN;
# }
warn "From:[$in] to [$out].*.gz\nFQin:[$fq1],[$fq2]\n";

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


sub getFQitem($) {
	my $FH = $_[0];
	my @dat;
	return undef if eof($FH);
	for (1 .. 4) {
		my $line = <$FH>;
		die '[x]FQ file not match 1 !' unless defined($line);
		push @dat,$line;
	}
	$dat[0] =~ /^\@([^\/]+)\/(\d+)$/ or die;
	my ($id,$r12) = ($1,$2);
	$id = (split /#/,$id)[0];
	return [$id,$r12,\@dat];
}

my $OUT1 = openpipe('gzip -9c',"$out.1.fq.gz");
my $OUT2 = openpipe('gzip -9c',"$out.2.fq.gz");

my ($READLEN,$InsSum,$InsCnt,%IDs,%InsDat)=(0,0,0);

sub Sam2FQ($) {
	my $rd = $_[0];
	my ($rd12);
	if ($$rd[1] & 0x40) {
		$rd12 = 1;
	} elsif ($$rd[1] & 0x80) {
		$rd12 = 2;
	}
	die "CIAGR:$$rd[5]" if $$rd[5] =~ /H/;
	my $seq = $$rd[9];
	my $qual = $$rd[10];
	if ($$rd[1] & 0x10) {
		$seq =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
		$seq = reverse $seq;
		$qual = reverse $qual;
	}
	my $id = join('','@',$$rd[0],'/',$rd12);
	my $str = join("\n",$id,$seq,'+',$qual);
	return [$rd12,$str];
}
sub doSamPair($$) {
	my ($rd1,$rd2) = @_;
	if ($inType eq 'BSMAP') {
		++$IDs{$$rd1[0]};
	} elsif ($inType eq 'bwa-meth') {
		my $ret1 = Sam2FQ($rd1);
		my $ret2 = Sam2FQ($rd2);
		($ret1,$ret2) = sort {$a->[0] <=> $b->[0]} ($ret1,$ret2);
		print $OUT1 $ret1->[1],"\n";
		print $OUT2 $ret2->[1],"\n";
	}
}
open( IN,"-|","samtools view -F768 $in") or die "Error opening $in: $!\n";
while (my $line = <IN>) {
	#my ($id, $flag, $ref, $pos, $mapq, $CIAGR, $mref, $mpos, $isize, $seq, $qual, @OPT) = split /\t/,$line;
	#print "$id, $flag, $ref, $pos, $mapq, $CIAGR, $mref, $mpos, $isize\n";
	my @Dat1 = split /\t/,$line;
	my $line2 = <IN>;
	die '[x]SAM/BAM file not paired !' unless defined($line2);
	my @Dat2 = split /\t/,$line2;
	#if ( $ref eq 'chrEBV' or ($flag & 12) ) {
	if ( $Dat1[2] eq 'chrEBV' or ($Dat1[1] & 12) or $Dat2[2] eq 'chrEBV' or ($Dat2[1] & 12) ) {
		if ($Dat1[0] ne $Dat2[0]) {
			die "[x]SAM/BAM file not paired as $Dat1[0] != $Dat2[0] !\n";
		}
		#++$IDs{$Dat1[0]};
		doSamPair(\@Dat1,\@Dat2);
		print $OUT0 "$line$line2";
	} elsif ("$Dat1[5]$Dat2[5]" =~ /\d+S/) {
		my $flag = 0;
		my (@edit_cmd) = "$Dat1[5]$Dat2[5]" =~ m/\d+S/g;
		foreach my $cmd (@edit_cmd) {
			if (my ($S) = $cmd =~ m/(\d+)S/) {
				if ($S >= $minSoftClipLen) {
					$flag = 1;
					last;
				}
			}
		}
		if ($flag == 1) {
			doSamPair(\@Dat1,\@Dat2);
			print $OUT0 "$line$line2";
		} else {
			print $OUT4 "$line$line2";
		}
	} elsif ( $Dat1[2] eq '*' or $Dat2[2] eq '*' ) {
		#++$IDs{$Dat1[0]};
		doSamPair(\@Dat1,\@Dat2);
		warn "-1->[$line$line2";
	} else {
		my $isize = abs($Dat1[8]);
		if ($Dat1[8] + $Dat2[8] == 0) {
			++$InsDat{$isize};
			if ($Dat1[8] != 0) {
				$InsSum += $isize;
				++$InsCnt;
			}
			if ($Dat1[5] =~ /^(\d+)M$/) {	# should be outside this `if`, but inside means less times thus faster(?).
				$READLEN = $1 if $READLEN < $1;
			}
		} else {
			if ($Dat1[8] == $Dat2[8]) {
				warn "-2->$isize\->[$line$line2";
			} else {
				warn "-3->$isize\->[$line$line2";
			}
		}
	}
}
close IN;
$InsCnt = -1 unless $InsCnt;
print OUT3 "BamIn:$in\nFQ1in:$fq1\nFQ2in:$fq2\nOut:$out\nReadLength:$READLEN\nInsertSize:",$InsSum/$InsCnt,"\n\n_Total_\t$InsCnt\n";
for my $s (sort {$a<=>$b} keys %InsDat) {
	print OUT3 "$s\t$InsDat{$s}\n";
}
close OUT3;
close $OUT0;
close $OUT4;

if ($inType eq 'BSMAP') {
	my $FQ1 = openfile($fq1);
	my $FQ2 = openfile($fq2);
	while (my $FQ1dat = getFQitem($FQ1)) {
		my $FQ2dat = getFQitem($FQ2);
		unless ( $FQ1dat->[0] eq $FQ2dat->[0] ) {
			die '[x]FQ file not match 2 ! ',$FQ1dat->[0],'|',$FQ2dat->[0];
		}
		if (exists $IDs{$FQ1dat->[0]}) {
			print $OUT1 join('',@{$FQ1dat->[2]});
			print $OUT2 join('',@{$FQ2dat->[2]});
			delete $IDs{$FQ1dat->[0]};
		}
		last unless keys %IDs;
		#warn "$FQ1dat->[0]\n";
	}
	close $FQ1;
	close $FQ2;
}
close $OUT1;
close $OUT2;

#system("samtools view -bS $out.sam.gz >$out.bam");
__END__
bsmap -u -z 64 -p 12 -v 10 -q 2 -d HomoGRCh38.fa -a F12HPCCCSZ0010_Upload/s00_C.bs_1.fq.gz -b F12HPCCCSZ0010_Upload/s00_C.bs_2.fq.gz -o s00_C.bs.bam 2> s00_C.bs.err
bsmap -u -z 64 -p 12 -v 10 -q 2 -d HomoGRCh38/HomoGRCh38.fa.gz -a F12HPCCCSZ0010_Upload/s01_P.bs_1.fq.gz -b F12HPCCCSZ0010_Upload/s01_P.bs_2.fq.gz -o s01_P.bs.bam 2> s01_P.bs.err
-rw-r--r-- 1 huxs users  8598234466 Dec 22 21:53 s00_C.bs.bam
-rw-r--r-- 1 huxs users         961 Dec 22 21:53 s00_C.bs.err
-rw-r--r-- 1 huxs users           0 Dec 22 16:32 s00_C.bs.log

./getUnPaired.pl s00_C.bs.bam n_grep
../getUnPaired.pl s00_C.bshum.bam
