#!/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

my (%gen_code,%t,@t,%start_codes,$code,%AA,@AAs);
$code=<<Ecode;
  AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
Starts = ---M---------------M---------------M----------------------------
Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
Ecode

@t=split /\n/,$code;
map {s/\s//g;@_=split /=/;$t{$_[0]}=[split //,$_[1]]} @t;
for (@{$t{AAs}}) {
	++$AA{$_};
	my $base=(shift @{$t{Base1}}).(shift @{$t{Base2}}).(shift @{$t{Base3}});
#print "$base";
	my $start=shift @{$t{Starts}};
	++$start_codes{$base} if $start eq 'M';
	$gen_code{$base}=$_;
#print " -> [$_]\n";
}
@AAs = sort keys %AA;

ddx \%gen_code;
ddx \%start_codes;
ddx scalar @AAs,\@AAs;

my @Bases = sort qw( A T C G );
my %mutBases;
$mutBases{$_} = [] for @Bases;
for my $mut (@Bases) {
	for (@Bases) {
		push @{$mutBases{$_}},$mut if $mut ne $_;
	}
}
ddx \%mutBases;

for my $pos (0..2) {
	my ($cnsSyn,$cntAAc) = (0,0);
	for my $codon (sort keys %gen_code) {
		my $ntApos = substr $codon,$pos,1;
		my $oriAA = $gen_code{$codon};
		next if $oriAA eq '*';
		print "$pos\t$codon:$oriAA\t";
		for my $mut ( @{$mutBases{$ntApos}} ) {
			my $mutCodon = $codon;
			substr $mutCodon,$pos,1,$mut;
			my $newAA = $gen_code{$mutCodon};
			next if $newAA eq '*';
			print "$mutCodon:$newAA";
			if ( $oriAA eq $newAA ) {
				++$cnsSyn;
				print ".";
			} else {
				++$cntAAc;
				print "x";
			}
			print " ";
		}
		print "\n";
	}
	print ">>$pos: $cnsSyn,$cntAAc,",$cnsSyn/($cnsSyn+$cntAAc),"\n";
}

__END__
$ perl get1986property.pl |grep \>\>
>>0: 8,166,0.0459770114942529
>>1: 0,176,0
>>2: 126,50,0.715909090909091


