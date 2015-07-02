#!/usr/bin/env perl
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );

my $SAMTOOLSBIN="samtools";

die "Usage: $0 <output> <sam/bam files>\n" if @ARGV <2;

my $out = shift;

warn "[!]Stat [$out] <- [",join('],[',@ARGV),"]\n";
open OUT,'>',$out or die "Error opening $out: $!\n";

my %StatDat;
for (@ARGV) {
	print STDERR "[!]Reading [$_]: ";
	open INSAM,'-|',"$SAMTOOLSBIN view $_" or die "Error opening $_: $!\n";
	close INSAM;
	warn "done.\n";
}

__END__
./bam.pl stat.log  /share/users/huxs/work/tiger2/part1/PD2M-LSJ-BHX01*
