#!/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120330
=cut
use strict;
use warnings;
#use Data::Dump qw(ddx);

die "Usage: $0 <blast annot file>\n" if @ARGV<1;
my $annot=shift;

my (%Stat,%StatF);
open H,'<',$annot or die "Error opening $annot : $!\n";
while (<H>) {
	next if /^#/;
	chomp;
	my ($qseqid,$sacc,$length,$evalue,$mismatch,$annot)=split /\t/;
	$annot =~ s/^\[//;
	$annot =~ s/\]$//;
#print "$qseqid,$sacc,$length,$evalue,$mismatch,$annot\n";
	my $isPREDICTED=0;
	if ($annot =~ s/^PREDICTED: //) {
		$isPREDICTED=1;
	}
	$qseqid=join(' ',(split / /,$annot)[0..1]);
	$sacc=(split /,/,$annot)[0];
#print "$annot -> $qseqid | $sacc\n";
	if (exists $Stat{$qseqid}) {
		++$Stat{$qseqid}->[0];
		$Stat{$qseqid}->[1] += $isPREDICTED;
	} else {
		$Stat{$qseqid} = [1,$isPREDICTED];
	}
	if (exists $StatF{$sacc}) {
		++$StatF{$sacc}->[0];
		$StatF{$sacc}->[1] += $isPREDICTED;
	} else {
		$StatF{$sacc} = [1,$isPREDICTED];
	}
}
close H;

open O,'>',$annot.'.stat' or die "Error opening $annot.stat : $!\n";
print O join("\t",split(' ',q/#Sp All_Count PREDICTE_Count/)),"\n";
for my $k (sort {$Stat{$b}->[0] <=> $Stat{$a}->[0] || $Stat{$a}->[1] <=> $Stat{$b}->[1] || $a cmp $b} keys %Stat) {
	print O join("\t",$k,@{$Stat{$k}}),"\n";
}
print O "\n#NEXT",'-' x 72,"\n\n";
for my $k (sort {$StatF{$b}->[0] <=> $StatF{$a}->[0] || $StatF{$a}->[1] <=> $StatF{$b}->[1] || $a cmp $b} keys %StatF) {
	print O join("\t",$k,@{$StatF{$k}}),"\n";
}
close O;
