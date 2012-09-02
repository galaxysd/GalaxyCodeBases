#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Data::Dump qw(ddx);
use Statistics::LineFit;
use Galaxy::IO;
use Galaxy::SeqTools;

die "Usage: $0 <in> <out.prefix>\n" if @ARGV<2;
my $inf=shift;
my $outf=shift;

my (%dat,$t,%out);
open I,'<',$inf or die $!;
while (<I>) {
	chomp;
	my @items = split /\s+/;
	# 0NMid 1chr 2+- 3start 4end 5gene 6incmpl5 7incmpl3 8scaffold 9start 10end 11+- 12gene(same)
	if ( $items[2] eq $items[11] ) {
		$t = 0;
	} else {
		$t = 1;
	}
	push @{$dat{$items[8]}{$items[1]}}, [$t,$items[5],$items[3],$items[4],$items[9],$items[10]];
}
close I;
#ddx \%dat;

for my $scafd (keys %dat) {
	my (%scoreSame,%scoreAbs,%thePos);
	for my $chr (keys %{$dat{$scafd}}) {
		my $alignmtsA = $dat{$scafd}{$chr};
		my (%cntStrand);
		@$alignmtsA = sort {$a->[2] <=> $b->[2]} @$alignmtsA;
		for my $item (@$alignmtsA) {
print join("\t",$scafd,$chr,@$item);
			++$cntStrand{$$item[0]};
			my $lenChr = $$item[3] - $$item[2];
			my $lenScafd = $$item[5] - $$item[4];
			my ($a,$b) = sort {$a<=>$b} ($lenChr,$lenScafd);
			if ($a) {
				$t = $a / $b;
			} else {
#print "\t***\n";
				#next;
				$t = 0.1;
			}
			push @{$thePos{$chr}},[$$item[0],$$item[2]+$lenChr/2,$$item[4]+$lenScafd/2];
print join("\t",'',$lenScafd,$lenChr,$t),"\n";
			if ($$item[0]) {
				$scoreSame{$chr} -= $t;
			} else {
				$scoreSame{$chr} += $t;
			}
			$scoreAbs{$chr} += $t;
		}
		print join(",",$scoreAbs{$chr},$scoreSame{$chr},%cntStrand),"\t-----\n";
	}
	my ($theChr) = sort { $scoreAbs{$b} <=> $scoreAbs{$a} } keys %scoreAbs;
	my $lineFit = Statistics::LineFit->new();
	my $diffStrand;
	if ($scoreSame{$theChr} > 0) {
		$diffStrand = 0;
	} else {
		$diffStrand = 1;
	}
	print $theChr,'-' x 70,"\n";
}

open O,'>',$outf.'.plst' or die $!;

close O;

__END__

