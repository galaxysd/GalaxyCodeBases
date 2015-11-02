#!/usr/bin/perl -w
use strict;
use warnings;
use Data::Dump;

open I,'<','resLOC.anno' or die $!;
open A,'<','resLOC.cydat' or die $!;
open O,'>','resLOC.ricexpro' or die $!;
open ON,'>','resLOC.rpratio' or die $!;

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
close I;

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
#   44020 => ["LOC_Os11g24450", "LOC_Os11g25030"],
#   44023 => ["LOC_Os05g37830"],

#ddx \%Loc2Acc;
#   LOC_Os12g37720 => { 32729 => 1, 34753 => 1, 41627 => 1 },
#   LOC_Os12g37770 => { 19813 => 1 },

#ddx \%Acc2Accession;	# undef exists from file !

my (%cyDat,%cyRDat);

my %RXP2Name = (
RXP_1006 => 'ABA', #Shoot gene expression profile in response to abscisic acid
RXP_1009 => 'BRs', #Shoot gene expression profile in response to brassinosteroid
RXP_1010 => 'CK' #Shoot gene expression profile in response to cytokinin
RXP_1008 => 'Aux' #Shoot gene expression profile in response to auxin
RXP_1007 => 'GA' #Shoot gene expression profile in response to gibberellin
RXP_1012 => 'JA' #Shoot gene expression profile in response to jasmonic acid
);	# Cy3 (mock treatment) and Cy5 (hormone treatment) => Cy5/Cy3, http://ricexpro.dna.affrc.go.jp/RXP_1006/index.php

my (@FeatureOrder, @ExpOrder, @FeatNameOrder, @reps, %rep, @FinalOrder, @RatioExpOrder, @RatioFinalOrder);
# FeatureNum    RXP_1006        RXP_1009        RXP_1010
## 0 hr (Cy3)|0 hr (Cy5)|1 hr (Cy3)|1 hr (Cy5)|3 hr (Cy3)|3 hr (Cy5)|6 hr (Cy3)|6 hr (Cy5)|12 hr (Cy3)|12 hr (Cy5)
while (<A>) {
	if ( /^#/ ) {
		#print O $_;
		if ( /^# FeatureNum\b/ ) {
			chomp;
			@FeatureOrder = split /\t/;
			shift @FeatureOrder;
			for (@FeatureOrder) {
				my $tmp = $RXP2Name{$_};
				push @FeatNameOrder,$tmp;
				print "$_, $tmp\n";
			}
		} elsif ( /^## (.+)$/ ) {
			$_ = $1;
			chomp;
			my @Exps = split /\|/;
			for ( @Exps ) {
				/^(\d+) hr \(Cy([35])\)$/ or die;
				push @ExpOrder,"${1}h_Cy$2";
				push @RatioExpOrder,"${1}h" if $2 eq '3';
				#print "$_ => ${1}h_Cy$2\n";
			}
			print "@ExpOrder -> @RatioExpOrder\n";
		}
		next;
	}
	my ( $acc,@dat ) = split /\t/;
	my (%DAT,%RDAT);
	for my $FNO ( @FeatNameOrder ) {
		my $tmp = shift @dat;
		my @rep = split /\,/,$tmp;
		for my $repline ( @rep ) {
			my @repdat = split /\|/,$repline;
			my $repid = shift @repdat;
			if ($repid eq 'mean') {
				$RDAT{$FNO} = \@repdat;
			} else {
				++$rep{$repid};
				$DAT{$FNO}{$repid} = \@repdat;
			}
		}
	}
	$cyDat{$acc} = \%DAT;
	$cyRDat{$acc} = \%RDAT;
}
close A;
@reps = sort keys %rep;
for my $FNO ( @FeatNameOrder ) {
	for my $RO ( @reps ) {
		for my $EO ( @ExpOrder ) {
			my $t = $RO;
			$t =~ s/^rep/r/;
			push @FinalOrder,"$FNO$EO$t";
		}
	}
	for my $REO ( @RatioExpOrder ) {
		push @RatioFinalOrder,"$FNO$REO";
	}
}
for my $acc ( keys %cyDat ) {
	my %cyhash = %{$cyDat{$acc}};
	my %cyrhash = %{$cyRDat{$acc}};
	my (@item,@ritem);
	for my $FNO ( @FeatNameOrder ) {
		for my $RO ( @reps ) {
			for my $EO ( 0 .. $#ExpOrder ) {
				push @item,$cyhash{$FNO}{$RO}->[$EO];
			}
		}
		for my $EO ( 0 .. $#RatioExpOrder ) {
			push @ritem,$cyrhash{$FNO}->[1+2*$EO] / $cyrhash{$FNO}->[2*$EO];
			#push @ritem,log($cyrhash{$FNO}->[1+2*$EO] / $cyrhash{$FNO}->[2*$EO])/log(2);
		}
	}
	$cyDat{$acc} = \@item;
	$cyRDat{$acc} = \@ritem;
}
#ddx \%cyDat;
print "\n@FinalOrder\n";
print join(',',@{$cyDat{'45209'}}),"\n";
print "\n@RatioFinalOrder\n";
print join(',',@{$cyRDat{'45209'}}),"\n";
#   45209 => [
#              "rep1|8.0073|3.6185|8.9315|5.0357|4.6752|4.7962|83.208|12.445|44.395|6.3024,rep2|3.1927|2.6697|3.3677|4.4205|3.3884|3.1185|7.425|6.9132|5.2888|5.7744,mean|5.6|3.1441|6.1496|4.7281|4.0318|3.9574|45.317|9.6791|24.842|6.0384",
#              "rep1|8.0073|3.6185|4.2184|3.7634|3.4171|6.2807|3.9199|4.3256|5.1421|3.9737,rep2|3.1927|2.6697|4.048|3.8231|4.0791|2.9637|3.569|4.0095|2.597|2.6025,mean|5.6|3.1441|4.1332|3.7933|3.7481|4.6222|3.7445|4.1676|3.8695|3.2881",
#              "rep1|8.0073|3.6185|3.2538|4.5685|11.197|3.8527|4.4721|3.4333|4.8344|4.8462,rep2|3.1927|2.6697|2.7079|3.9205|7.4918|3.104|3.1729|2.6352|6.9551|4.1439,mean|5.6|3.1441|2.9809|4.2445|9.3444|3.4784|3.8225|3.0343|5.8947|4.495\n",
#            ],
#   45209 => {
#              ABA => {
#                       rep1 => [
#                                 8.0073,
#                                 3.6185,
#                                 8.9315,
#                                 5.0357,
#                                 4.6752,
#                                 4.7962,
#                                 83.208,
#                                 12.445,
#                                 44.395,
#                                 6.3024,
#                               ],
#                       rep2 => [
#                                 3.1927,
#                                 2.6697,
#                                 3.3677,
#                                 4.4205,
#                                 3.3884,
#                                 3.1185,
#                                 7.425,
#                                 6.9132,
#                                 5.2888,
#                                 5.7744,
#                               ],
#                     },
#              BRs => {
#                       rep1 => [
#                                 8.0073,
#                                 3.6185,
#                                 4.2184,
#                                 3.7634,
#                                 3.4171,
#                                 6.2807,
#                                 3.9199,
#                                 4.3256,
#                                 5.1421,
#                                 3.9737,
#                               ],
#                       rep2 => [
#                                 3.1927,
#                                 2.6697,
#                                 4.048,
#                                 3.8231,
#                                 4.0791,
#                                 2.9637,
#                                 3.569,
#                                 4.0095,
#                                 2.597,
#                                 2.6025,
#                               ],
#                     },
#              CK  => {
#                       rep1 => [
#                                 8.0073,
#                                 3.6185,
#                                 3.2538,
#                                 4.5685,
#                                 11.197,
#                                 3.8527,
#                                 4.4721,
#                                 3.4333,
#                                 4.8344,
#                                 4.8462,
#                               ],
#                       rep2 => [
#                                 3.1927,
#                                 2.6697,
#                                 2.7079,
#                                 3.9205,
#                                 7.4918,
#                                 3.104,
#                                 3.1729,
#                                 2.6352,
#                                 6.9551,
#                                 4.1439,
#                               ],
#                     },
#            },
# 8.0073,3.6185,8.9315,5.0357,4.6752,4.7962,83.208,12.445,44.395,6.3024,3.1927,2.6697,3.3677,4.4205,3.3884,3.1185,7.425,6.9132,5.2888,5.7744,8.0073,3.6185,4.2184,3.7634,3.4171,6.2807,3.9199,4.3256,5.1421,3.9737,3.1927,2.6697,4.048,3.8231,4.0791,2.9637,3.569,4.0095,2.597,2.6025,8.0073,3.6185,3.2538,4.5685,11.197,3.8527,4.4721,3.4333,4.8344,4.8462,3.1927,2.6697,2.7079,3.9205,7.4918,3.104,3.1729,2.6352,6.9551,4.1439

print O join("\t",'YORF','NAME',@FinalOrder),"\n";
for my $loc (sort keys %Loc2Acc) {
	for my $acc ( sort keys %{$Loc2Acc{$loc}} ) {
		print O join("\t",$loc,$acc,@{$cyDat{$acc}}),"\n";
	}
}
close O;
print ON join("\t",'YORF','NAME',@RatioFinalOrder),"\n";
for my $loc (sort keys %Loc2Acc) {
	for my $acc ( sort keys %{$Loc2Acc{$loc}} ) {
		print ON join("\t",$loc,$acc,@{$cyRDat{$acc}}),"\n";
	}
}
close ON;
__END__
