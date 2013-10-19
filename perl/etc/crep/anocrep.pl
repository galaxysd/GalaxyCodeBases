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

my %outLOC;
open I,'<',$inf or die;
while (<I>) {
	chomp;
	my @dat = split /\t/;
	my $aid = $dat[0];
	if (exists $AFF2LOC{$aid}) {
		for my $i ( @{ $AFF2LOC{$aid} } ) {
			$outLOC{ $i->[0] } = [0,[]] unless exists $outLOC{ $i->[0] };
			push @{ $outLOC{$i->[0]}->[1] },[@dat];
			#push @{ $outLOC{$i->[0]}->[1] },\@dat;
			$outLOC{$i->[0]}->[0] += $dat[1] * $i->[2];
		}
	}
}
close I;
#ddx \%outLOC;
#   LOC_Os01g02300 => [
#                       2.99555110455193,
#                       [
#                         [
#                           "Os.27046.1.S1_at",
#                           2.48051948051948,
#                           101.87,
#                           41.07,
#                           "66|189.6|50",
#                           "49.5|59.1|14.6",
#                         ],
#                         [
#                           "Os.27046.2.S1_x_at",
#                           3.46376167185416,
#                           259.67,
#                           74.97,
#                           "88.3|399.4|291.3",
#                           "13.9|167|44",
#                         ],
#                       ],
#                     ],
#   LOC_Os01g02300 => {
#                       "Os.27046.1.S1_at"   => [10, 0.476190476190476],
#                       "Os.27046.2.S1_x_at" => [11, 0.523809523809524],
#                     },

open O,'>',$inf.'.out' or die;
for my $tid (sort { $outLOC{$b}->[0] <=> $outLOC{$a}->[0] } keys %outLOC) {
	my @items = @{ $outLOC{$tid} };
	my @affprobes;
	for ( @{ $items[1] } ) {
		$$_[1] = int(0.5+1000*$$_[1])/1000;
		push @affprobes,join('~',@{$LOC2AFF{$tid}{$$_[0]}},@$_);
	}
	print O join("\t",$tid,(scalar @affprobes),$items[0],$LOCdesc{$tid},join("\t",@affprobes) ),"\n";
}
close O;

__END__
perl anocrep.pl crep_anno_merged.tsv crep_all_tsv_new.txt.up2
perl anocrep.pl crep_anno_merged.tsv crep_all_tsv_new.txt.down2
