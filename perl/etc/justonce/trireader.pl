#!/usr/bin/perl
use strict;
use warnings;
use List::Util 1.26 qw(sum0);
use Data::Dump qw(ddx);

my $fhead = 'snp/bam.head';
my @files = qw[snp/SNP.out snp/tomor1.snp.out snp/tumour2.snp.out];

sub getDat($) {
	my $fh = $_[0];
	while(<$fh>) {
		chomp;
		my @dat = split /\t/;
		next unless $dat[6] eq 'PASS';
		my @dep1 = split /,/,$dat[9];
		my @dep2 = split /,/,$dat[10];
		my $sum1 = sum0(@dep1);
		my $sum2 = sum0(@dep2);
		next unless ($sum1+$sum2) > 20;
		#return [@dat[0,1,7,6,9,10],$sum1+$sum2];
		return [@dat[0,1,7]];
	}
	return ["\t",-1,"NN"];
}

my (@ChrIDs,%ChrLen,%Chrank);
open H,'<',$fhead or die $!;
while (<H>) {
	chomp;
	my @dat = split /\t/;
	next if $dat[0] ne '@SQ';
	my $cid = (split /:/,$dat[1])[1];
	my $clen = (split /:/,$dat[2])[1];
	push @ChrIDs,$cid;
	$ChrLen{$cid} = $clen;
}
close H;
%Chrank = map { $ChrIDs[$_], $_ + 1 } 0 .. $#ChrIDs;
#ddx \%Chrank,\@ChrIDs,\%ChrLen;

my ($fha,$fhb,$fhc);
open $fha,'<',$files[0] or die $!;
open $fhb,'<',$files[1] or die $!;
open $fhc,'<',$files[2] or die $!;
<$fha>;<$fha>;<$fhb>;<$fhb>;<$fhc>;<$fhc>;

my @IDs = qw[A B C];
my %FHs = (
	'A' => [$fha,getDat($fha)],
	'B' => [$fhb,getDat($fhb)],
	'C' => [$fhc,getDat($fhc)],
);
#print $FHs{'A'}->[1][0],"--\n";

open O,'>','res.tsv' or die $!;
while (($FHs{'A'}->[1][0] ne "\t") and ($FHs{'B'}->[1][0] ne "\t") and ($FHs{'C'}->[1][0] ne "\t")) {
	#my %ChrIDs = map { $_ => $FHs{$_}->[1][0] } @IDs;
	my @aID = sort { $Chrank{$FHs{$b}->[1][0]} <=> $Chrank{$FHs{$a}->[1][0]} || $FHs{$b}->[1][1] <=> $FHs{$a}->[1][1] } @IDs;	# desc
	#ddx \%FHs,\@aID;
	for my $i (1 .. $#aID) {
		while (
			(($Chrank{$FHs{$aID[0]}->[1][0]} == $Chrank{$FHs{$aID[$i]}->[1][0]}) and ($FHs{$aID[0]}->[1][1] > $FHs{$aID[$i]}->[1][1]))
			or ($Chrank{$FHs{$aID[0]}->[1][0]} > $Chrank{$FHs{$aID[$i]}->[1][0]})
		) {
			$FHs{$aID[$i]}->[1] = getDat($FHs{$aID[$i]}->[0]);
		}
	}
	if ($FHs{'A'}->[1][1] == $FHs{'B'}->[1][1] and $FHs{'A'}->[1][1] == $FHs{'C'}->[1][1]) {
		#ddx \%FHs;
		my $type='N';
		if ( $FHs{'A'}->[1][2] eq $FHs{'B'}->[1][2] and $FHs{'A'}->[1][2] eq $FHs{'C'}->[1][2] ) {
			$type='AAA';
		} elsif ( $FHs{'A'}->[1][2] eq $FHs{'B'}->[1][2] and $FHs{'A'}->[1][2] ne $FHs{'C'}->[1][2] ) {
			$type='AAB';
		} elsif ( $FHs{'A'}->[1][2] ne $FHs{'B'}->[1][2] and $FHs{'A'}->[1][2] eq $FHs{'C'}->[1][2] ) {
			$type='ABA';
		} elsif ( $FHs{'A'}->[1][2] ne $FHs{'B'}->[1][2] and $FHs{'B'}->[1][2] eq $FHs{'C'}->[1][2] ) {
			$type='ABB';
		} elsif ( $FHs{'A'}->[1][2] ne $FHs{'B'}->[1][2] and $FHs{'A'}->[1][2] ne $FHs{'C'}->[1][2] and $FHs{'B'}->[1][2] ne $FHs{'C'}->[1][2] ) {
			$type='ABC';
		}
		print O join("\t",$FHs{'A'}->[1][0],$FHs{'A'}->[1][1],$type,join(',',$FHs{'A'}->[1][2],$FHs{'B'}->[1][2],$FHs{'C'}->[1][2])),"\n";
		for my $f (@aID) {
			$FHs{$f}->[1] = getDat($FHs{$f}->[0]);
		}
	} else {
		$FHs{'A'}->[1] = getDat($FHs{'A'}->[0]);
	}
}
close O;
close $fha;close $fhb;close $fhc;
