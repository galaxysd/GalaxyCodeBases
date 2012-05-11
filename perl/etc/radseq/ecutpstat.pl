#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120330
=cut
use strict;
use warnings;
#use Data::Dump qw(ddx);

die "Usage: $0 <reference> <out>\n" if @ARGV<2;
my $fa=shift;
my $out=shift;

my $Eseq="CTGCAG";
my $EcutAt=5;

open O,'>',$out or die "Error opening $out : $!\n";
print O "# Ref: [$fa]\n#Enzyme: [$Eseq], Cut after $EcutAt\n\n#ChrID\tCut\n";
my ($CountAll,%Count);

open I,'<',$fa or die "Error opening $fa : $!\n";
while (<I>) {
    s/^>//;
	/^(\S+)/ or next;
	my $seqname = $1;
    #print STDERR " >$seqname ...";
	$/=">";
	my $genome=<I>;
	chomp $genome;
	$genome=~s/\s//g;
	$/="\n";
	#print STDERR "\b\b\b   \t",length $genome,".\n";

	while ( $genome =~ m/$Eseq/g ) {
		# push @ret, [ $-[0], $+[0] ];
		print O join("\t",$seqname,$-[0]+$EcutAt),"\n";	# starts from 1, so OK to use directly for "cut after".
		++$CountAll;
		++$Count{$seqname};
	}

	$genome='';
}
close I;

print O "\nCut stat:\n";
for (sort keys %Count) {
	print O "# $_\t$Count{$_}\n"
}
print O "\n#Total: $CountAll";
close O;
