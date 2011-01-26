#!/bin/env perl
use strict;
use warnings;
use lib '/nas/RD_09C/resequencing/soft/lib';
#use GalaxyXS::ChromByte;
use Data::Dump qw(ddx);

unless (@ARGV > 0) {
    print "perl $0 <markerpos file>\n";
    exit 0;
}

my ($file)=@ARGV;
open I,'<',$file or die "Error:[$file] $!\n";
open O,'>',$file.'.f' or die "Error:[$file.f] $!\n";
$_=<I>;
print O $_;
while (<I>) {
	next if /^#/;
	chomp;
	my ($Qid,$mcM,$Sid,$pos,$strand,$Pidentity,$E,$BTOP)=split /\t/;
	next if $E > 3e-11;
	print O join("\t",$Qid,$mcM,$Sid,$pos,$strand,$Pidentity,$E,$BTOP),"\n";
}
close I;
close O;

__END__
find ./markerpos/ -name '*.pos'|xargs -n1 ./markerposfilter.pl
