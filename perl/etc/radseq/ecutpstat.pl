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
my ($CountAll,%Count,%Len);

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

	$Count{$seqname} = 0;
	$Len{$seqname} = length $genome;
	while ( $genome =~ m/$Eseq/g ) {
		# push @ret, [ $-[0], $+[0] ];
		print O join("\t",$seqname,$-[0]+$EcutAt),"\n";	# starts from 0, so OK to use directly with $EcutAt for "cut after".
		++$CountAll;
		++$Count{$seqname};
	}

	$genome='';
}
close I;

my @notCut;
print O "\n#Cut stat:\n#ChrID\tLen\tCut\tAvg_distance\n";
for (sort keys %Count) {
	print O '# ',join("\t",$_,$Len{$_},$Count{$_},$Count{$_}?int(0.5 + $Len{$_}/(1+$Count{$_})):'.'),"\n";
	push @notCut,$_ unless $Count{$_};
}
print O "#Total: $CountAll\n#notCut: ",scalar @notCut,"\n";
print O "#notCut:[",join(',',@notCut),"]\n" if @notCut;
close O;
