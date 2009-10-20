#!/bin/env perl
use strict;
use warnings;

# 20090805 init
# 20090807 patch for NULL sequence output. Any reads shorter than 30
#  will lengthern by 10 with tailing 'N' x 10 .
# 20090807 patch again for a minimal output length of 30.
# 20091019 print stat info to STDERR

# for 1M reads, 280M vf required.
unless (@ARGV){
	print "perl $0 <adapter list file> <fq file>\n";
	exit;
}

my $adapter = shift;
my $fq = shift;

my %adapter;
open LIST, "$adapter" || die "$!\n";
while (<LIST>) {
	chomp;
	if (/^FC/) {
		my @line = split /\s+/;
		if (defined $adapter{$line[0]}) {
			unless ($adapter{$line[0]}>$line[2]) {
				$adapter{$line[0]} = $line[2];
			}
		} else {
			$adapter{$line[0]} = $line[2];
		}
	}
}
close LIST;

my ($maxreadlen,$filtedreads,$copyedreads,$inbp,$outbp,$readlen)=(0);

my $read_num = 0;
open FQ, "$fq" || die "$!\n";
while (my $line1=<FQ>) {
	chomp $line1;
	my $line2 = <FQ>;
	chomp $line2;
	$readlen=length $line2;
	$maxreadlen = $readlen if $maxreadlen < $readlen;
	$inbp += $readlen;
	my $line3 = <FQ>;
	chomp $line3;
	my $line4 = <FQ>;
	chomp $line4;
	my ($read_name) = $line1=~ /^\@(FC.*)/;
	++$read_num;
	if (defined $adapter{$read_name}) {
		++$filtedreads;
#		$line1 =~ s/FC\w{9}:\d+:\d+:\d+:\d+/_${read_num}_/;
		my $line2_new = substr($line2, 0, $adapter{$read_name});
		my $line4_new = substr($line4, 0, $adapter{$read_name});
		my $ii=30-length($line2_new);
		$outbp += length($line2_new);	# well, every perl variable is a object, thus length takes no much time.
		if ($adapter{$read_name} < 30) {
			$line2_new .= 'N' x $ii;
			$line4_new .= 'B' x $ii;
		}
		print "$line1\n$line2_new\n$line3\n$line4_new\n";
	} else {
		++$copyedreads;
		$outbp += $readlen;
	#	$line1 =~ s/FC\w{9}:\d+:\d+:\d+:\d+/_${read_num}_/;
		print "$line1\n$line2\n$line3\n$line4\n";
	}
}
close FQ;

warn "# $read_num parsed in [$fq]
# Using [$adapter]
# MaxReadLen\t$maxreadlen
# InReads\t$read_num
# InBPs\t$inbp
# FiltedReads\t$filtedreads
# CopyedReads\t$copyedreads
# OutBP\t$outbp
# All done !\n";
