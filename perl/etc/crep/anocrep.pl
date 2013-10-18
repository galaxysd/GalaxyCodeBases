#!/bin/env perl
use strict;
use warnings;
use IO::Unread qw(unread);
use Data::Dump qw(ddx);

die "Usage: $0 <anno_file> <input>\n" if @ARGV < 2;
my ($anno,$inf)=@ARGV;

my (%ANNO,%ANNOcnt,%LOCdesc,%LOC2AFF,%AFF2LOC);

open I,'<',$anno or die;
while (<I>) {
	chomp;
	my @dat = split /\t/;
	next unless exists $dat[1];
	next if $dat[1] eq 'NONE|NONE|NONE';
	my $id = shift @dat;
	my $sum=0;
	my @tids;
	for (@dat) {
		my ($tid,$ratio,$desc) = split /\|/;
		$ratio =~ s#/11$## or die "[$ratio]";
		$sum += $ratio;
		if (exists $LOCdesc{$tid} and $LOCdesc{$tid} ne $desc) {
			die "[x]$_\n$LOCdesc{$tid}\n";
		}
		$LOCdesc{$tid} = $desc;
		push @tids,[$tid,$ratio];
	}
	$ANNO{$id} = \@tids;
	$ANNOcnt{$id} = $sum;
}
close I;
#ddx \%ANNO;
#   "Os.10031.1.S1_at"          => [["LOC_Os05g35140", 11], ["LOC_Os07g12230", 1]],

for my $aid (keys %ANNO) {
	my @items = @{$ANNO{$aid}};
	for (@items) {
		$LOC2AFF{$_->[0]} = [[],0] unless exists $LOC2AFF{$_->[0]};
=pod
		my @dat = @{ $LOC2AFF{$_->[0]} };
		push @{ $dat[0] },$_->[0];
		$dat[1] += $_->[1];
		$LOC2AFF{$_->[0]} = \@dat;
=cut
		push @{ $LOC2AFF{$_->[0]}->[0] },[$aid,$_->[1]];
		$LOC2AFF{$_->[0]}->[1] += $_->[1];
	}
}
#ddx \%LOC2AFF;
#   LOC_Os05g35140 => [[["Os.10031.1.S1_at", 11]], 11],
#   LOC_Os07g12230 => [
#                       [
#                         ["OsAffx.22438.1.S1_x_at", 11],
#                         ["Os.10031.1.S1_at", 1],
#                         ["OsAffx.22438.1.S1_at", 11],
#                         ["OsAffx.8676.1.S1_at", 4],
#                       ],
#                       27,
#                     ],

for my $tid (keys %LOC2AFF) {
	my @dat = @{$LOC2AFF{$tid}->[0]};
	my $sum = $LOC2AFF{$tid}->[1];
	my %hdat;
	for (@dat) {
		$_->[2] = $_->[1] / $sum;
		$hdat{ $_->[0] } = [$_->[1],$_->[2]];
		$AFF2LOC{ $_->[0] } = [] unless exists $AFF2LOC{ $_->[0] };
		push @{ $AFF2LOC{$_->[0]} },[$tid,$_->[1],$_->[2]];
	}
	$LOC2AFF{$tid} = \%hdat;
}
#ddx \%LOC2AFF;
#   LOC_Os05g35140 => { "Os.10031.1.S1_at" => [11, 1] },
#   LOC_Os07g12230 => {
#                       "Os.10031.1.S1_at"       => [1, 0.037037037037037],
#                       "OsAffx.22438.1.S1_at"   => [11, 0.407407407407407],
#                       "OsAffx.22438.1.S1_x_at" => [11, 0.407407407407407],
#                       "OsAffx.8676.1.S1_at"    => [4, 0.148148148148148],
#                     },
#ddx \%AFF2LOC;
#   "Os.10031.1.S1_at"          => [
#                                    ["LOC_Os05g35140", 11, 1],
#                                    ["LOC_Os07g12230", 1, 0.037037037037037],
#                                  ],

die;
my %outLOC;
open I,'<',$inf or die;
while (<I>) {
	chomp;
	my @dat = split /\t/;
	my $aid = shift @dat;
	if (exists $AFF2LOC{$aid}) {
		for my $i ( @{ $AFF2LOC{$aid} } ) {
			$outLOC{ $i->[0] } = [] unless exists $outLOC{ $i->[0] };
		}
	}
}
close I;

__END__
perl anocrep.pl crep_anno_merged.tsv crep_all_tsv_new.txt.up2
