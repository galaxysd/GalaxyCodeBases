#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Data::Dump qw(ddx);
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
	push @{$dat{$items[1]}{$items[8]}}, [$t,$items[5],$items[3],$items[4],$items[9],$items[10]];
}
close I;
#ddx \%dat;

for my $chr (keys %dat) {
	for my $scafd (keys %{$dat{$chr}}) {
		my $alignmtsA = $dat{$chr}{$scafd};
		my ($scoreSame,%cntStrand)=(0);
		for my $item (@$alignmtsA) {
print join("\t",$chr,$scafd,@$item);
			++$cntStrand{$$item[0]};
			my $lenChr = $$item[3] - $$item[2];
			my $lenScafd = $$item[5] - $$item[4];
			my ($a,$b) = sort {$a<=>$b} ($lenChr,$lenScafd);
			if ($a) {
				$t = $a / $b;
			} else {
				next;
				$t = 0;
			}
print join("\t",'',$lenChr,$lenScafd,$t),"\n";
			if ($$item[0]) {
				$scoreSame += $t;
			} else {
				$scoreSame -= $t;
			}
		}
print join(",",$scoreSame,%cntStrand),"\n",'-'x75,"\n";
	}
}

open O,'>',$outf.'.plst' or die $!;

close O;

__END__

