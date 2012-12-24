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
my (%Scaff,@Scaffs);
while (<L>) {
	chomp;
	push @Scaffs,$_ unless exists $Scaff{$_};
	++$Scaff{$_};
}
close L;
print "Total ",scalar keys %Scaff," Scaffolds: [",join('],[',sort keys %Scaff),"]\n";

open LOG,'>',"$outfs.log" or die "Error opening $outfs.log : $!\n";
print LOG "Total ",scalar keys %Scaff," Scaffolds: [",join('],[',sort keys %Scaff),"]\n";

$t=1;
for (@Scaffs) {
	$Scaff{$_} = $t;
	print LOG "$_: $t\n";
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

for my $type (qw[case control all]) {
	my $inprefix;
	if ($type eq 'all') {
		$inprefix = $infs;
	} else {
		$inprefix = "${infs}.$type";
	}
	for my $scaffid (keys %Scaff) {
		open C,'<',"$inprefix.tped" or die "Error opening $inprefix.tped : $!\n";
		open O,'>',"$outfs.$type.$scaffid.tped" or die "Error opening $outfs.$type.$scaffid.tped : $!\n";
		#symlink "$inprefix.tfam","$outfs.$type.$scaffid.tfam";
		while (<C>) {
			my @items = split /\t/;
			if (exists $Rid{$items[1]} and $Rid{$items[1]} eq $scaffid) {
				$items[0] = $Scaff{$Rid{$items[1]}};
				print O join("\t",@items);
			}
		}
		close C;
		close O;
		my $cmd = "p-link --tped \'$outfs.$type.$scaffid.tped\' --tfam $inprefix.tfam --out \'${outfs}$type.$scaffid\' --recode";
		print "Running [$cmd].\n";
		print LOG "$cmd\n";
		system($cmd);
	}
}
close LOG;

__END__

./ped2hap.pl tighap.lst radseq17 dohap17

p-link --tfile dohap17.case --out dohap17case --recode
