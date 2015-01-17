#!/usr/bin/perl -w
#
#
use strict;

my $flen = shift || 500;
my $rlen = shift || 80;
my $n_rd = shift || 20;
my $m_ol = shift || 20;
my $m = shift || 1000;

for(my $i=0;$i<$m;$i++){
	my @ss = ();
	for(my $j=0;$j<$n_rd;$j++){
		push(@ss, int($flen * rand()));
	}
	@ss = sort {$a <=> $b} @ss;
	my $len = $rlen;
	for(my $j=0;$j<$n_rd;$j++){
		last if($ss[$j] + $m_ol > $len);
		$len = $ss[$j] + $rlen;
	}
	print "$len\n";
}

1;
