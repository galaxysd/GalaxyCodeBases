#!/usr/bin/perl -w
use strict;
use warnings;
use Data::Dump;

#die "Usage: $0 <s> <t>\nW(AA), W(Aa) and W(aa) is 1-s, 1, 1-t respectively.\n" if @ARGV<2;

my $totalLen = 1000;

my $mfa=<<Emfa;
1  TAGACAATGACATACGTTATGTTGAA
2  ..T...TA.G.C.C...C..T.....
3  ..T...T.C..C.C.AA..CT.....
4  .TTCACT...GCG.A...T.T.G.T.
5  C.T...T..G.C.C...C..TC.C.G
Emfa
my @t=split /\n/,$mfa;
my %dat;
map {@_=split /\s+/;$dat{$_[0]}=[split //,$_[1]]} @t;
#ddx \%dat;
my @seqid = sort keys %dat;
my $seqcnt = @seqid;
my $seqlen = @{$dat{'1'}};

print "$_: ",join('',@{$dat{$_}}),"\n" for @seqid;
print "=== @seqid ===\n$mfa";

my @PolyT;
for my $i (0 .. $#{$dat{'1'}}) {
	my $ref = $dat{'1'}->[$i];
	for my $k (keys %dat) {
		next if $k eq '1';
		if ($dat{$k}->[$i] eq '.') {
			$dat{$k}->[$i] = $ref;
			#$PolyT[$i] = 0;
		} else {
			$PolyT[$i] = 1;
		}
	}
}
print "P: ",join('',@PolyT),"\n\n";
my $PolyC = 0;
for (@PolyT) {
	++$PolyC if $_ and $_ eq '1';
}

print "$_: ",join('',@{$dat{$_}}),"\n" for @seqid;
print "=== $seqcnt,$seqlen,$PolyC ===\n";

sub getTpi(@) {
	print "--I: @_   ";
	my %t;
	++$t{$_} for @_;
	my @t = sort {$a <=> $b} values %t;
	die if @t != 2;
	my ($p,$q) = @t;
	my $sum = $p + $q;
	$p /= $sum; $q /= $sum;
	my $Tpi = 2*$p*$q*$seqcnt / ($seqcnt-1);
	print "$p,$q,$Tpi\n";
	return $Tpi;
}

my $Tpi = 0;
for my $i (0 .. $#{$dat{'1'}}) {
	my @arr;
	for my $k (@seqid) {
		push @arr,$dat{$k}->[$i];
	}
	$Tpi += getTpi(@arr);
}
print $Tpi;
$Tpi /= $totalLen;
print ", $Tpi\n";

my $wUnder = 0;
for my $i ( 1 .. ($seqcnt-1) ) {
	$wUnder += 1/$i;
	print "$i: $wUnder\n";
}
my $Tw = $PolyC / $wUnder / $totalLen;
print "-> $Tw, $wUnder\n";

print "Tpi: $Tpi   Tw: $Tw\n";


