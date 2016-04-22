#!/usr/bin/perl
use strict;
use warnings;

open I, "-|", "find ../b.bwa/ -name *.sam.gz";
open O, ">", "makefile";

my %sam;
while (<I>) {
	chomp;
	/^\.\.\/b\.bwa\/(\w{10})\/(\w{12})\.sam\.gz$/;
	$sam{$1}{$2} = 1;
}
close I;

print O "bam:";
print O " $_.bam" foreach sort keys %sam;
print O "\n";

foreach my $a (sort keys %sam) {
	$a =~ /(^\w{6})/;
	my $sm = $1;
	my $ss;
	foreach my $b (sort keys %{$sam{$a}}) {
		my $g = "../b.bwa/${a}/$b.sam.gz";
		my $u = "${a}_$b.bam";
		my $s = "${a}_$b.sort";
		$ss .= " $s.bam";
		print O "$s.bam: $u\n\tsamtools sort $u $s\n\trm -f $u\n";
		print O "$u: $g\n\t./addRG.pl \"\@RG ID:$b PL:ILLUMINA SM:$sm\" $g $u\n";
	}
	print O "$a.bam:$ss\n\tsamtools merge -$ss |samtools rmdup - $a.bam 2>$a.rmdup.log\n\trm -f $ss\n";
}
close O;

