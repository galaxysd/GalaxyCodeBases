#!/usr/bin/env perl
use strict;
use warnings;

my ($in,$out) = @ARGV;

sub revcom($) {
	my $str = $_[0];
	$str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
	my $rev = reverse $str;
	$rev =~ tr/[](){}<>/][)(}{></;
	return $rev;
}

sub openfile($) {
	my ($filename)=@_;
	my $infile;
	if ($filename=~/.xz$/i) {
		open( $infile,"-|","xz -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.gz$/i) {
		open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.bz2$/i) {
		open( $infile,"-|","bzip2 -dc $filename") or die "Error opening $filename: $!\n";
	} else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
	return $infile;
}

my $infh = openfile($in);
my $outfh;
if ($in=~/.gz$/i) {
	open( $outfh,"|-","gzip -c $out") or die "Error opening $out: $!\n";
} else {
	open $outfh,'>',$out or die "Error opening $out: $!\n";
}
