#!/bin/env perl
use strict;
use warnings;

# 20090805 init
# 20090807 patch for NULL sequence output. Any reads shorter than 30
#  will lengthern by 10 with tailing 'N' x 10 .
# 20090807 patch again for a minimal output length of 30.
# 20091019 print stat info to STDERR

# for 1M reads, 280M vf required.
# See also: /nas/RD_09C/resequencing/soft/pipeline/index_pipe/bin/filter_fq_peZ.pl

unless (@ARGV){
	print "perl $0 <adapter list file> <fq file> <out.fq.gz> <nfo> <max len>\n";
	exit;
}

my ($adapter,$fq,$out,$nfo,$maxL)=@ARGV;

my %adapter;
#open LIST,'<',$adapter || die "$!\n";
if ($adapter =~ /\.gz$/) {
	open( LIST,'-|',"gzip -dc $adapter") or die "[x]Error on $adapter: $!\n";
} elsif ($adapter =~ /\.bz2$/) {
 	open( LIST,'-|',"bzip2 -dc $adapter") or die "[x]Error on $adapter: $!\n";
} else {open( LIST,'<',$adapter) or die "[x]Error on $adapter: $!\n";}
while (<LIST>) {
	next if /^#/;
	chomp;
	my @line = split /\s+/;
	if (defined $adapter{$line[0]}) {
		$adapter{$line[0]} = $line[2] if $adapter{$line[0]} < $line[2];
	} else {
		$adapter{$line[0]} = $line[2];
	}
}
close LIST;

my ($maxreadlen,$filtedreads,$copyedreads,$inbp,$outbp,$readlen)=(0);

my $read_num = 0;
#open FQ, "$fq" || die "$!\n";
if ($fq =~ /\.gz$/) {
	open( FQ,'-|',"gzip -dc $fq") or die "[x]Error on $fq: $!\n";
} elsif ($fq =~ /\.bz2$/) {
 	open( FQ,'-|',"bzip2 -dc $fq") or die "[x]Error on $fq: $!\n";
} else {open( FQ,'<',$fq) or die "[x]Error on $fq: $!\n";}

while (my $line1=<FQ>) {
	chomp $line1;
	my $line2 = <FQ>;
	chomp $line2;
	$readlen=length $line2;
#	$maxreadlen = $readlen if $maxreadlen < $readlen;
	$inbp += $readlen;
	my $line3 = <FQ>;
	chomp $line3;
	my $line4 = <FQ>;
	chomp $line4;
	my ($read_name) = $line1=~ /^\@(.*)$/;	# /^\@(\S+)/
	++$read_num;
	if (defined $adapter{$read_name}) {
		++$filtedreads;
#		$line1 =~ s/\@\S+\#/\@\_${read_num}_\#/;
		my $alen=$adapter{$read_name};
		#my $line2_new = substr($line2, 0, $adapter{$read_name});
		substr($line2, $alen, $readlen-$alen, 'N' x ($readlen-$alen) );
		substr($line4, $alen, $readlen-$alen, 'B' x ($readlen-$alen) );
	#	$outbp += $alen;
		$readlen = $alen;
		#my $ii=30-length($line2_new);
		#$outbp += length($line2_new);	# well, every perl variable is a object, thus length takes no much time.
		#if ($adapter{$read_name} < 30) {
		#	$line2_new .= 'N' x $ii;
		#	$line4_new .= 'B' x $ii;
		#}
	#	print "$line1\n$line2\n$line3\n$line4\n";
	} else {
		++$copyedreads;
	#	$outbp += $readlen;
	#	$line1 =~ s/FC\w{9}:\d+:\d+:\d+:\d+/_${read_num}_/;
	#	print "$line1\n$line2\n$line3\n$line4\n";
	}

	if ($maxL && $maxL<$readlen) {
		$line2=substr $line2,0,$maxL;
		$line4=substr $line4,0,$maxL;
		$readlen=$maxL;
	}
	$maxreadlen = $readlen if $maxreadlen < $readlen;

	print "$line1\n$line2\n$line3\n$line4\n";
	$outbp += $readlen;
}
close FQ;

warn "$read_num parsed in [$fq]
Using [$adapter]
MaxReadLen\t$maxreadlen
InReads\t$read_num
InBPs\t$inbp
FiltedReads\t$filtedreads
CopyedReads\t$copyedreads
OutReads\t",$filtedreads+$copyedreads,"
OutBP\t$outbp
All done !\n";
