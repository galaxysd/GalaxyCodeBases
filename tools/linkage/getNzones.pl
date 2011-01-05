#!/bin/env perl
use strict;
use warnings;

unless (@ARGV > 0) {
    print "perl $0 <fa file> <out N zones>\n";
    exit 0;
}

my ($in,$out)=@ARGV;

my $infile;
if ($in =~ /.bz2$/) {
	open( $infile,"-|","bzip2 -dc $in") or die "Error: $!\n";
} elsif ($in =~ /.gz$/) {
 	open( $infile,"-|","gzip -dc $in") or die "Error: $!\n";
} else {open( $infile,"<",$in) or die "Error: $!\n";}

{
	local $/=">";
	$_=<$infile>;
	die "[x]Not a FASTA file !\n" unless /^\s*>/;
}

open O,'>',$out or die "Error: $!\n";
while (<$infile>) {
	chomp;
	my ($id,$desc)=split / /,$_,2;
	$desc='' unless $desc;
	$/=">";
	my $seq=<$infile>;
	chomp $seq;
	$seq =~ s/\s//g;
	$/="\n";
	print STDERR ">$id,\t[$desc]:\n";
	my $len=length($seq);
	my $pos=0;
	my $inNzone=0;
	my ($Nstart,$Ncount,$Nlen)=(0,0,0);
	while ($pos <= $len) {
		my $base=substr($seq,$pos,1);
		my $isN = ($base =~ /N/i);
		if (! $isN and ! $inNzone) {
			$inNzone=0;
		} elsif ($isN and ! $inNzone) {
			$inNzone=1;
			$Nstart=$pos;
			print O "$id\t",$pos+1,"\t";
		} elsif (! $isN and $inNzone) {
			$inNzone=0;
			print O join("\t",$pos,$pos-$Nstart),"\n";
			$Nlen += $pos-$Nstart;
			++$Ncount;
			$Nstart=0;
		} else {
			;
		}
		++$pos;
	}
	my $avg='NA';
	$avg=int($Nlen/$Ncount) if $Ncount;
	print STDERR " Len:$len\tN-zone: $Nlen <- $Ncount x ",$avg,"\n";
}
close O;
close $infile;
