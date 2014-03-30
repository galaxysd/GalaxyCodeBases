#!/usr/bin/perl -w
use strict;
use warnings;
use Data::Dump;

#die "Usage: $0 <s> <t>\nW(AA), W(Aa) and W(aa) is 1-s, 1, 1-t respectively.\n" if @ARGV<2;

my @N = (1000,1000,100,1000,1000);
ddx \@N;

sub getNe($$) {
	my ($m,$f) = @_;
	my $Ne = 4*$m*$f/($m+$f);
	return $Ne;
}

sub getPop($$) {
	my ($Pm,$Pf) = @_;
	my $mf = $Pm + $Pf;
	my $sum = 0;
	print STDERR '--- Ne: ';
	for my $n ( @N ) {
		my $Ne = getNe( $n*$Pm/$mf,$n*$Pf/$mf );
		$sum += 1/$Ne;
		print STDERR "$Ne, ";
	}
	print STDERR "\n";
	$sum /= scalar @N;
	return 1/$sum;
}

print STDERR "1:1\n",getPop(1,1),"\n";

print STDERR "1:4\n",getPop(1,4),"\n";
