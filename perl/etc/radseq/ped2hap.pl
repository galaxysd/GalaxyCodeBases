#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
Purpose: Read bcf, get tped for p-link
Notes: rad2marker is deprecated.
=cut
use strict;
use warnings;
use Data::Dump qw(ddx);
use Data::Dumper;

die "Usage: $0 <hap scaffold list> <tped sets prefix> <out>\n" if @ARGV<3;
my $lstfs=shift;
my $infs=shift;
my $outfs=shift;

my (%Stat,$t);
open L,'<',$lstfs or die "Error opening $lstfs : $!\n";
my %Scaff;
while (<L>) {
	chomp;
	++$Scaff{$_};
}
close L;
print "Total ",scalar keys %Scaff," Scaffolds: [",join('],[',sort keys %Scaff),"]\n";
$t=1;
for (sort keys %Scaff) {
	$Scaff{$_} = $t;
	++$t;
}
die if $t > 22;

open D,'<',$infs.'.dict' or die "Error opening $infs.dict : $!\n";
my %Rid;
while (<D>) {
	chomp;
	my ($scaffid,$pos,$rid) = split /\t/;
	if (exists $Scaff{$scaffid}) {
		$Rid{$rid} = $scaffid;
	}
}
close D;
print "Total ",scalar keys %Rid," Markers.\n";

for my $type (qw[case control]) {
	open C,'<',"$infs.$type.tped" or die "Error opening $infs.$type.tped : $!\n";
	open O,'>',"$outfs.$type.tped" or die "Error opening $outfs.$type.tped : $!\n";
	symlink "$infs.$type.tfam","$outfs.$type.tfam";
	while (<C>) {
		my @items = split /\t/;
		if (exists $Rid{$items[1]}) {
			$items[0] = $Scaff{$Rid{$items[1]}};
			print O join("\t",@items);
		}
	}
	close C;
	close O;
}

__END__

./ped2hap.pl tighap.lst radseq17 dohap17
