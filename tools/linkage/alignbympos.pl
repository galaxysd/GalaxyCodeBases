#!/bin/env perl
use strict;
use warnings;
use lib '/nas/RD_09C/resequencing/soft/lib';
#use GalaxyXS::ChromByte;
use Data::Dump qw(ddx);

unless (@ARGV > 0) {
    print "perl $0 <markerpos scaffold> <markerpos chr> <out>\n";
    exit 0;
}

my ($scaff,$chrf,$outf)=@ARGV;
open S,'<',$scaff or die "Error:[$scaff] $!\n";

open C,'<',$chrf or die "Error:[$chrf] $!\n";
