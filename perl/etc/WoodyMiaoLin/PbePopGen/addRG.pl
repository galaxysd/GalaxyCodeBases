#!/usr/bin/perl
use strict;
use warnings;

die "This program add HEADER and TAG of READ GROUP to bam file.\nAuther: Woody\nUsage: $0 <in.bam> <out.bam>\n" if @ARGV<2;

my $in = shift;
my $out = shift;
$in =~ /^.+\/(\w+)\.rmdup\.bam$/;
my $ip = $1;
$out =~ /^(.+)\.bam$/;
my $op = $1;

open I, "-|", "samtools-0.1.7 view -h $in";
open HEAD, ">", "$op.header.sam";
open READ, "|-", "samtools-0.1.7 view -bS - >$op.read.bam";

warn "$ip adding TAG...\n";
my %RG;
while (<I>) {
	if (/^@/) {
		print HEAD;
		print READ;
	} else {
		chomp;
		my @a = split /\t/;
		my @b = split /:/, $a[0];
		my $rg = "$b[0]_$b[1]_$b[2]_$b[3]";
		$RG{$rg} = 1;
		print READ "$_\tRG:Z:$rg\n";
	}
}

close I;
close READ;

foreach (sort keys %RG) {
	print HEAD "\@RG\tID:$_\tPL:ILLUMINA\tSM:$ip\n";
}
close HEAD;

warn "$ip reheadering...\n";
system "samtools reheader -P $op.header.sam $op.read.bam >$op.bam";
system "rm $op.header.sam $op.read.bam";
warn "Done!\n";
