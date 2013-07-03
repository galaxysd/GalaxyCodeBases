#!/bin/env perl
use strict;
use warnings;

die "Usage: $0 <col> <input[gzipped]> <output(gzipped)>\n" if @ARGV < 3;
my ($col,$inf,$outf)=@ARGV;
--$col;

my ($Sum,$Cnt)=(0,0);

sub openFH($) {
	my $inf = $_[0];
	my $FH;
	if ($inf =~ /\.gz$/i) {
		open $FH,'-|',"gzip -dc $inf" or die "Error opening $inf: $!\n";
	} elsif ($inf =~ /\.xz$/i) {
		open $FH,'-|',"xz -dc $inf" or die "Error opening $inf: $!\n";
	} elsif ($inf =~ /\.bz2$/i) {
		open $FH,'-|',"bzip2 -dc $inf" or die "Error opening $inf: $!\n";
	} else {
		open $FH,'<',$inf or die "Error opening $inf: $!\n";
	}
	return $FH;
}

my $FH = openFH($inf);
while (<$FH>) {
	my @Dat = split /\s+/;
	next if /^[#;]/;
	$Sum += $Dat[$col];
	++$Cnt;
}
close $FH;
print "Input $Cnt lines, Sum = $Sum\n";

$FH = openFH($inf);
open OUT,'|-', "gzip -9c >$outf" or die "Error opening $outf: $!\n";
while (<$FH>) {
	my @Dat = split /\s+/;
	if (/^[#;]/) {
		print OUT $_;
		next;
	}
	print OUT join("\t",@Dat,$Dat[$col]/$Sum),"\n";
}
close $FH;
close OUT;
