#!/bin/env perl
use strict;
use warnings;

die "Usage: $0 <art.aln file>" unless @ARGV == 1;

open I,'<',$ARGV[0] or die;
my ($ref,$seq,@seq,@ref,%Mis);
my ($BaseCount,$MisBaseCount) = (0,0);
while(<I>) {
	chomp($ref=<I>);
	chomp($seq=<I>);
	#print "[$ref] [$seq]\n";
	@seq=split //,$seq;
	@ref=split //,$ref;
	my $mis=0;
	for (@seq) {
		my $base=shift @ref;
		++$mis if $base ne $_;
	}
	++$Mis{int(0.5+100*$mis/(length $seq))};
	$BaseCount += length $seq;
	$MisBaseCount += $mis;
	#print "[$ref] [$seq] $mis\n";
}

print "mis%\tcount\n",'-'x75,"\n";
for (sort {$a<=>$b} keys %Mis) {
	print "$_\t$Mis{$_}\n";
}
print "\nMisRate: $MisBaseCount/$BaseCount = ",$MisBaseCount/$BaseCount,"\n";
