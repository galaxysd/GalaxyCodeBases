#!/bin/env perl
use strict;
use warnings;

# 20090805 init
# 20090807 patch for NULL sequence output. Any reads shorter than 30
#  will lengthern by 10 with tailing 'N' x 10 .
# 20090807 patch again for a minimal output length of 30.

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

my $read_num = 0;
open FQ, "$fq" || die "$!\n";
while (my $line1=<FQ>) {
	chomp $line1;
	my $line2 = <FQ>;
	chomp $line2;
	my $line3 = <FQ>;
	chomp $line3;
	my $line4 = <FQ>;
	chomp $line4;
	my ($read_name) = $line1=~ /^\@(FC.*)/;
	if (defined $adapter{$read_name}) {
		++$read_num;
#		$line1 =~ s/FC\w{9}:\d+:\d+:\d+:\d+/_${read_num}_/;
		my $line2_new = substr($line2, 0, $adapter{$read_name});
		my $line4_new = substr($line4, 0, $adapter{$read_name});
		my $ii=30-length($line2_new);
		if ($adapter{$read_name} < 30) {
			$line2_new .= 'N' x $ii;
			$line4_new .= 'B' x $ii;
		}
		print "$line1\n$line2_new\n$line3\n$line4_new\n";
	} else {
		++$read_num;
	#	$line1 =~ s/FC\w{9}:\d+:\d+:\d+:\d+/_${read_num}_/;
		print "$line1\n$line2\n$line3\n$line4\n";
	}
}
close FQ;

warn "$read_num parsed in [$fq]\nUsing [$adapter]\n";
