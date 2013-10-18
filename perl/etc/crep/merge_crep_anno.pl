#!/bin/env perl
use strict;
use warnings;
#use IO::Unread qw(unread);
use Term::ANSIColor qw(:constants);
use Data::Dump qw(ddx);

die "Usage: $0 <out> <input1> [input2 ...]\n" if @ARGV < 2;
my $out = shift @ARGV;

my $psif='./Rice.psi';

my (@psid,$tmp,%ANNO,%INPUTIDS);
open PSI,'<',$psif or die "Error: $!\n";
while (<PSI>) {
	next if /^#/;
	$tmp = (split /\t/,$_)[1];
	push @psid,$tmp;
	++$INPUTIDS{$tmp};
	#print "[$tmp_line]\n";
}
close PSI;
$tmp = $#psid + 1;
warn GREEN,BOLD,$tmp,RESET," ProbeSets loaded.\n";

while (<>) {
	chomp;
	my @dat = split /\t/;
	next unless exists $dat[1];
	$ANNO{shift @dat} = \@dat;
}
$tmp = (keys %ANNO);
warn GREEN,BOLD,$tmp,RESET," ProbeSets Found.\n";

if ($tmp < $#psid + 1) {
	print "Lacking:\n";
	for (@psid) {
		unless (exists $ANNO{$_}) {
			print "[$_]\n";
		}
	}
	die;
}
if ($tmp > $#psid + 1) {
	print "More:\n";
	for (keys %ANNO) {
		unless (exists $INPUTIDS{$_}) {
			print "[$_]\n";
		}
	}
}

open O,'>',$out or die "Error: $!\n";
for (@psid) {
	print O join("\t",$_,@{$ANNO{$_}}),"\n";
}
close O;

__END__
perl merge_crep_anno.pl crep_anno_merged.tsv crep_anno_tsv.all.txt crep_anno_tsv.txt*
wc -l crep_anno_merged.tsv
