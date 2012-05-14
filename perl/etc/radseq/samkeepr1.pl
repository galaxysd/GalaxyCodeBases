#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120330
=cut
use strict;
use warnings;
#use Data::Dump qw(ddx);

die "Usage: $0 <sam in> <out prefix>\n" if @ARGV<2;
my $in=shift;
my $outp=shift;

sub openfile($) {
    my ($filename)=@_;
    my $infile;
    if ($filename=~/.bz2$/) {
	    open( $infile,"-|","bzip2 -dc $filename") or die "Error opening $filename: $!\n";
    } elsif ($filename=~/.gz$/) {
     	open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
    } else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
    return $infile;
}

my $samin = openfile($in);
open O,'|-',"gzip -9c >$outp.sam.gz" or die "Error opening $outp.sam.gz with gzip: $!\n";
open L,'>',$outp.'.log' or die "Error opening $outp.log: $!\n";
select(L);
$|=1;
print L "From [$in] to [$outp.sam.gz]\n";
my ($Total,$Out,$notOut)=(0,0,0);

while (<$samin>) {
	if (/^@\w\w\t\w\w:/) {
		print O $_;
		next;
	}
	my @read1=split /\t/;
	if (($read1[1] & 64) == 64) {
		print O $_;
		++$Out;
	} elsif (($read1[1] & 128) == 128) {
		++$notOut;
	}
	++$Total;
}

close $samin;
close O;
print L "Read_1: $Out , ",$Out/$Total,"\nRead_2: $notOut , ",$notOut/$Total,"\nTotal: $Total\nRemain: ",$Total-$Out-$notOut,"\n";
close L;
