#!/bin/env perl
use strict;
use warnings;
use integer;

die "Usage: $0 <fq list> <outfile>\n" if @ARGV != 2;
my ($in,$out)=@ARGV;

open IN,'<',$in or die "Error opening $in: $!\n";
open OUT,'>',$out or die "Error opening $out: $!\n";

my $i=2;
while (<IN>) {
	next if /^$/;
	++$i;
	chomp;
	my $full = $_;
	s/\.(fastq|fq)(\.gz)?$// or die "[$_]";
	my $extra = $2;
	my @path = split /\//;
	my $name = pop @path;
#warn "$_\n[$name][@path]\n";
	my $isRADSEQ = 'N';
	$isRADSEQ = 'Y' if $name =~ /_NoIndex_/;
	my $id = (split /[_-]/,$name)[2];
	$name =~ /(_R|\.)([12])(_|.)([A-Za-z0-9])*?$/ or die "[$name]";
#warn "$2\n";
	my $tmpstr = join('_',$id,($i-$2)/2).".$2";
	#$full = `readlink -nf $full`;
	print OUT join("\t",$tmpstr,$isRADSEQ,$full),"\n";
	if ($isRADSEQ eq 'Y') {
		symlink $full,"work/radseq/$tmpstr.fq$extra";
	} else {
		symlink $full,"work/parents/$tmpstr.fq$extra";
	}
}
close IN;
close OUT;

