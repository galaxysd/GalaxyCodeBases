#!/usr/bin/perl -w
use strict;

my $ref = shift or die("Usage: $0 <ref> <anytag> <id>\n");
my $at  = shift or die("Usage: $0 <ref> <anytag> <id>\n");
my $id  = shift;
die("Usage: $0 <ref> <anytag> <id>\n") unless(defined $id);

open(IN, $at) or die;
open(OUT, ">lr.fa") or die;
open(ANY, ">lr") or die;
my $f = 0;
my $i = 0;
while(<IN>){
	if(/^T\t(\d+)/){
		last if($f);
		$f = ($1 == $id);
	}
	if($f){
		print ANY $_;
		my @ts = split;
		print OUT ">", ($i * 2), "\n$ts[3]\n";
		print OUT ">", ($i * 2 + 1), "\n$ts[4]\n";
		$i ++;
	}
}
close IN;
close OUT;
close ANY;

`blat -noHead $ref lr.fa lr.fa.psl`;
print `best_blat_hits.pl lr.fa.psl | msort -k 14,n16 | cut -f 14,16,17 | OVERLAP=10 simu_asm.pl`;

1;
