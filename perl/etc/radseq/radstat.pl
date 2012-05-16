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

my $Eseq="CTGCAG";
my $EcutAt=5;
my $EfwdTerm=1-$EcutAt;
my $ErevTerm=0;

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
	    open( $infile,"-|","samtools view -h -F 132 $filename") or die "Error opening $filename: $!\n";
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

my (%eCut,%eDat,@ChrOrder);
while (<$ecin>) {
	next if /^(#|$)/;
	my ($chr,$pos) = split /\t/;
	push @{$eCut{$chr}},$pos;
	push @ChrOrder,$chr unless exists $eCut{$chr};
	$eDat{$chr.'\t'.$pos}=[0,0,0,0,0,0,0,0];	# [SumF,CountF, SumR,CountR] for average depth
	# first 4, Unique(XT:A:U); second 4, Repeat(XT:A:R).
	# XT:A:N is fragmental like "61S12M2D3M2D10M14S".
	# Mate-sw(XT:A:M) is either short or with many(like 7) mismatches(including 'N' on reads)
}
close $ecin;

my ($Total,$Out,$notOut)=(0,0,0);

my (%ChrLen);
while (<$samin>) {
	if (/^@\w\w\t\w\w:/) {
		print O $_;
		if (/^\@SQ\tSN:(\S+)\tLN:(\d+)$/) {
			if (exists $eCut{$1}) {
				$ChrLen{$1} = $2;
				print STDERR "Chr:[$1]\tLen:[$2], Cut:[",scalar @{$eCut{$1}},"]\n";
			} else {
				warn "Chr:[$1], Len:[$2] not cut.\n";
			}
		}
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
