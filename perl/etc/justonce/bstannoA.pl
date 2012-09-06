#!/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120330
=cut
use strict;
use warnings;
#use Data::Dump qw(ddx);

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

$|=1;
open O,'>',$blastf.'.anot' or die "Error opening $blastf.anot : $!\n";
print O "# Count qseqid sacc length evalue mismatch annot bitscore qstart qend sstart send btop sseqid\n";

sub deal($) {
	my $Aref=$_[0];
	my (%Dat,@Arr);
	for (@$Aref) {
		push @{$Dat{$_->[2]}},$_;
	}
	for my $k (keys %Dat) {
		my $Emin=(sort {$b->[10] <=> $a->[10] || $a->[7] <=> $b->[7]} @{$Dat{$k}})[0];	# len Desc, E Asc
		push @Arr,[$Emin,scalar @{$Dat{$k}}];
	}
	for (@Arr) {
		my $Emin=$_->[0];
		my ($qseqid,$sseqid,$sacc,$qstart,$qend,$sstart,$send,$evalue,$bitscore,$score,$length,$pident,$nident,$mismatch,$positive,$gapopen,$gaps,$btop)=@$Emin;
		my $annot=$sseqid;
		$annot =~ s/^gi//;
		if (exists $Anno{$annot}) {
			$annot = '['.$Anno{$annot}.']';
			print O join("\t",$_->[1],$qseqid,$sacc,$length,$evalue,$mismatch,$annot,$bitscore,$qstart,$qend,$sstart,$send,$btop,$sseqid),"\n";
		} else {
			print STDERR '^';
warn "$annot\n";
		}
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

