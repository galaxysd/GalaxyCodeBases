#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ PKU <galaxy001@gmail.com>
Version: 1.0.0 @ 20120530
=cut
use strict;
use warnings;
#use Data::Dump qw(ddx);

die "Usage: $0 <sam1> <sam2> <out prefix>\n" if @ARGV<2;
my $in1=shift;
my $in2=shift;
my $outp=shift;

sub openfile($) {
    my ($filename)=@_;
    my ($infile,$lastline);
    if ($filename=~/.bz2$/) {
	    open( $infile,"-|","bzip2 -dc $filename") or die "Error opening $filename: $!\n";
    } elsif ($filename=~/.gz$/) {
     	open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
    } else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
	my $line;
	while (defined($line=<$infile>)) {
		if ($line =~ /^@\w\w\t\w\w:/) {
			next;
		}
		chomp $line;
		my @items = split /\t/,$line;
warn "[",join(" | ",@items),"]\n";
		if ($items[1] & 16) {	# reverse
			next;
		}
		return [$infile,\@items];
	}
	die "File Error.";
}

my $samin1 = openfile($in1);
my $samin2 = openfile($in2);
open O,'|-',"gzip -9c >$outp.cmpsam.gz" or die "Error opening $outp.cmpsam.gz with gzip: $!\n";
open L,'>',$outp.'.log' or die "Error opening $outp.log: $!\n";
select(L);
$|=1;
select(STDOUT);
print L "From [$in1,$in2] to [$outp.cmpsam.gz]\n";
my ($Total,$Out,$notOut)=(0,0,0);

sub getNextRead($) {
	my $FHa=$_[0];
	my $t=$FHa->[1];
	$FHa->[1] = undef;
	my $line;
	while (defined($line = readline($$FHa[0]))) {
		chomp $line;
		my @items = split /\t/,$line;
		if ($items[1] & 128) {	# second read
			next;
		}
		$FHa->[1] = \@items;
		return $t;
	}
}

my ($arr1,$arr2);
while (1) {
	$arr1 = getNextRead($samin1) or last;
	$arr2 = getNextRead($samin2) or last;
	die "$$arr1[0] ne $$arr2[0]" if $$arr1[0] ne $$arr2[0];
	print "[$#$arr1] [$#$arr2] $$arr1[0]\n";
}

close $samin1->[0];
close $samin2->[0];
close O;
print L "\n";
close L;
