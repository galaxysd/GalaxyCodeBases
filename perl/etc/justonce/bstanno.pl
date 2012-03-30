#!/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 1.0.0 @ 20110803
=cut
use strict;
use warnings;
use Data::Dump qw(ddx);

die "Usage: $0 <nt header> <blast file>\n" if @ARGV<2;
my $header=shift;
my $blastf=shift;

my %Anno;
open H,'<',$header or die "Error opening $header : $!\n";
while (<H>) {
	chomp;
	s/^>// or warn "[x][$_]\n";
	my @dat=split /\001/;
	for (@dat) {
		/^gi(\S*) (.*)/ or die;
		$Anno{$1}=$2;
#print "[@dat] ----- $1,$2\n";
	}
}
close H;
warn "[!]Load Hearers done.\n";

open O,'>',$blastf.'.anot' or die "Error opening $blastf.anot : $!\n";
print O "# qseqid sacc annot length evalue mismatch bitscore qstart qend sstart send btop sseqid\n";

sub deal($) {
	my $Aref=$_[0];
	my $Emin=(sort {$b->[10] <=> $a->[10] || $a->[7] <=> $b->[7]} @$Aref)[0];	# len Desc, E Asc
ddx $Aref;
print '-'x5,"\n";
ddx $Emin;
print '-'x75,"\n";
	my ($qseqid,$sseqid,$sacc,$qstart,$qend,$sstart,$send,$evalue,$bitscore,$score,$length,$pident,$nident,$mismatch,$positive,$gapopen,$gaps,$btop)=@$Emin;
	my $annot=$sseqid;
	$annot = s/^gi//;
	if (exists $Anno{$annot}) {
		$annot = '['.$Anno{$annot}.']';
		print O join("\t",$qseqid,$sacc,$annot,$length,$evalue,$mismatch,$bitscore,$qstart,$qend,$sstart,$send,$btop,$sseqid),"\n";
	} else {
		print STDERR '^';
	}
	return;
}

my $lastID='';
my @lastDat;
open B,'<',$blastf or die "Error opening $blastf : $!\n";
while (<B>) {
	chomp;
	my @t=split /\t/;
	if ($t[0] eq $lastID) {
		push @lastDat,\@t;
	} else {
		&deal(\@lastDat) if @lastDat>0;
		@lastDat = (\@t);
		$lastID = $t[0];
	}
}

close B;
close O;

