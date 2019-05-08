#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 1.0.0 @ 20180516
=cut
use strict;
use warnings;
use Parse::CSV;

use lib '.';

use Data::Dump qw(ddx);
#use FGI::GetCPI;

my $pI = '5';

die "Usage: $0 <info.csv> <fam.csv> <chip path> [out path]\n" if @ARGV<3;
my ($fninfo,$fnfam,$pchip,$pout) = @ARGV;
$pout = '.' unless $pout;
die "[x]Cannot read info.csv [$fninfo].\n" unless -r -s $fninfo;
die "[x]Cannot read fam.csv [$fnfam].\n" unless -r -s $fnfam;
die "[x]Cannot read Chip Path [$pchip].\n" unless -r $pchip;
die "[x]Cannot use out Path [$pout].\n" unless -r -w -x $pout;
warn "[1]Info:[$fninfo], Fam:[$fnfam], CHIP:[$pchip] to Out:[$pout]\n";

my $cinfo = Parse::CSV->new(file => $fninfo, names => 1);
my $cfam = Parse::CSV->new(file => $fnfam, names => 1);
mkdir "$pout/0lst";
open O,'>',"$pout/0lst/fq.lst";

my %fqIngo;
while ( my $value = $cinfo->fetch ) {
	next if $value->{Cell} eq '';
	$fqIngo{$value->{Sample}} = [$value->{Cell},$value->{Lane},$value->{Index}];
}
die $cinfo->errstr if $cinfo->errstr;
#ddx \%fqIngo;

for (sort keys %fqIngo) {
	my @d = @{$fqIngo{$_}};
	print O join("\t",$_,join('/',$d[0],$d[1],join('_',$d[0],$d[1],$pI.$d[2]))),"\n";
}
close O;

__END__
./genlst.pl info.csv fam.csv ./chip

while ( my $value = $cfam->fetch ) {
	next if $value->{Child} eq '';
	#ddx $value;
}
die $cfam->errstr if $cfam->errstr;
