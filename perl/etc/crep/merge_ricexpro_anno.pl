#!/usr/bin/perl -w
use strict;
use warnings;
use Data::Dump;

open I,'<','resLOC.anno' or die $!;
open A,'<','resLOC.cydat' or die $!;
open O,'>','resLOC.ricexpro' or die $!;

my (%Acc2Loc,%Loc2Acc,%Acc2Accession);
while (<I>) {
	my ($tigLOC,undef,$FeatureNums,$Accessions) = split /\t/;
	#next unless defined $FeatureNums;
	my @FeatureNums = split /\|/,$FeatureNums;
	next unless @FeatureNums > 0;
#	my @Accessions = split /\|/,$Accessions;
#print "$tigLOC,@FeatureNums,@Accessions\n";
	++$Acc2Loc{$_}{$tigLOC} for @FeatureNums;
	++$Loc2Acc{$tigLOC}{$_} for @FeatureNums;
#	for my $i (0 .. $#FeatureNums) {
#		$Acc2Accession{$FeatureNums[$i]} = $Accessions[$i];
#	}
}

my @AccNums = sort {$a <=> $b} keys %Acc2Loc;
for my $acc (@AccNums) {
	my @LOCs = sort keys %{$Acc2Loc{$acc}};
	for (@LOCs) { die if $Acc2Loc{$acc}{$_}>1; }
	$Acc2Loc{$acc} = \@LOCs;
}
my @LocNums = sort keys %Loc2Acc;
for my $loc (@LocNums) {
	my @ACCs = sort keys %{$Loc2Acc{$loc}};
	for (@ACCs) {
		die if $Loc2Acc{$loc}{$_}>1;
		$Loc2Acc{$loc}{$_} = 1 / scalar(@{$Acc2Loc{$_}});
	}
	#$Loc2Acc{$loc} = \@ACCs;
}

#ddx \%Acc2Loc;
#ddx \%Loc2Acc;
#ddx \%Acc2Accession;	# undef exists from file !


__END__
