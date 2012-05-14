#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120330
=cut
use strict;
use warnings;
#use Data::Dump qw(ddx);

die "Usage: $0 <Cut sites> <sorted bam> <out prefix>\n" if @ARGV<2;
my $inec=shift;
my $insam=shift;
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
sub opensam($) {
    my ($filename)=@_;
    my $infile;
    if ($filename=~/.bam$/) {
	    open( $infile,"-|","samtools view -h -F 128 $filename") or die "Error opening $filename: $!\n";
    } elsif ($filename=~/.sam.gz$/) {
     	open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
    } elsif ($filename=~/.sam$/) {
     	open( $infile,"<",$filename) or die "Error opening $filename: $!\n";
	} else {die "[x]Only support .sam(.gz) & .bam [$filename]\n";}
    return $infile;
}

my $ecin = openfile($inec);
my $samin = opensam($insam);
open O,'|-',"gzip -9c >$outp.edep.gz" or die "Error opening $outp.edep.gz with gzip: $!\n";
open L,'>',$outp.'.edep.log' or die "Error opening $outp.edep.log: $!\n";
select(L);
$|=1;
print L "From [$samin],[$ecin] to [$outp.edep.gz]\n";

my %eCut;
while (<$ecin>) {
	next if /^#/;
	my ($chr,$pos) = split /\t/;
	push @{$eCut{$chr}},$pos;
}
close $ecin;

my ($Total,$Out,$notOut)=(0,0,0);

my (@ChrID,%ChrLen);
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
print L "Read_1: $Out , ",$Out/$Total,"\nRead_2: $notOut , ",$notOut/$Total,"Total: $Total\nRemain: ",$Total-$Out-$notOut,"\n";
close L;
