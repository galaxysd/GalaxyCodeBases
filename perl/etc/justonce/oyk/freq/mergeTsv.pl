#!/usr/bin/env perl
use strict;
use warnings;

use Data::Dump qw(ddx);

die "Usage: $0 <in1.tsv> <in2.tsv> ...\n[!]For duplicated rsIDs, later file favoried.\n" if @ARGV<1;
my (%Markers,%MarkerAF,$sum);
while (<>) {
	chomp;
	my ($rs,$chr,$pos,@d) = split /\t/,$_;
	$MarkerAF{$rs} = {@d};
	#$sum = eval join '+',@AFs;
	#print "$rs: $sum\n" if $sum != 1;
	$Markers{$rs} = [$chr,$pos];
}
#ddx \%Markers;
my @ChrIDs = map {"chr$_"} (1..22,'X','Y','MT','M');
my $i = 0;
my %L = map { $_ => $i++ } @ChrIDs;
my @rsids = sort {
	no warnings 'uninitialized';
	exists $L{$Markers{$a}->[0]} <=> exists $L{$Markers{$b}->[0]} ||
	$L{$Markers{$a}->[0]} <=> $L{$Markers{$b}->[0]} ||
	$Markers{$a}->[1] <=> $Markers{$b}->[1] ||
	$a cmp $b
} keys %Markers;
for my $id (@rsids) {
	my @d = @{$Markers{$id}};
	my @allels = sort { $MarkerAF{$id}{$a} <=> $MarkerAF{$id}{$b} || $a cmp $b } keys %{$MarkerAF{$id}};
	my @af = map {"$_\t$MarkerAF{$id}{$_}"} @allels;
	print join("\t",$id,@d,@af),"\n";
}
__END__
./mergeTsv.pl MGIEasy.freq.tsv NS18.gt.dp30.freq.tsv >nippt.both.tsv
