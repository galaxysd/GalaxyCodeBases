#!/usr/bin/perl
use strict;
use warnings;

die "\nAdd HEADER and TAG of READ GROUP to bam file.\nAuther: Woody\nUsage: $0 <in.bam> <out.bam>\n\n" if @ARGV<2;

my $in = shift;
my $out = shift;
$in =~ /^(.+).bam$/;
my $ip = $1;
$out =~ /^(.+).bam$/;
my $op = $1;

open I, "-|", "samtools view -h $in";
open HEAD, ">", "$op.header.sam";
open READ, "|-", "samtools view -bS - >$op.reads.bam";

warn "Adding TAG...\n";
my %RG;
while (<I>) {
	if (/^@/) {
		print HEAD;
		print READ;
	} else {
		chomp;
		my @a = split /\t/;
		my @b = split /:/, $a[0];
		my $rg = "${ip}_$b[2]_$b[3]";
		$RG{$rg} = 1;
		print READ "$_\tRG:Z:$rg\n";
	}
}

close I;
close READ;

foreach (sort keys %RG) {
	print HEAD "\@RG\tID:$_\tPL:ILLUMINA\tSM:$ip\n";
}

warn "Printing HEADER...\n";
my $time = localtime();
print HEAD "\@CO\tHEADER and TAG of READ GROUP added at $time. ID: FilePre_FlowId_LaneNum, SM: FilePre.\n";
close HEAD;

warn "Reheadering...\n";
system "samtools reheader $op.header.sam $op.reads.bam >$op.bam";

warn "Removing $op.header.sam $op.reads.bam...\n";
system "rm $op.header.sam $op.reads.bam";
warn "Done!\n";
