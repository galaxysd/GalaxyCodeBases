#!/bin/env perl
use strict;
use warnings;
use lib '/nas/RD_09C/resequencing/soft/lib';
#use GalaxyXS::ChromByte;
use Data::Dump qw(ddx);

unless (@ARGV > 0) {
    print "perl $0 <group_file> <tmap output> <Out_prefix>\n";
    exit 0;
}

my ($in,$mapf)=@ARGV;
my (%MarkerGroup,@MarkerOrder,%MarkerCM);

open I,'<',$in or die $!;
{
	local $/=">\n";
	while (<I>) {
		chomp;
		my @dat=split /\n/;
		#print join("|",@dat),"\n";
		my ($id,$name)=split /=/,shift @dat,2;
		#print "$id\t$name\n";
		die "[x]ID confict found !" if exists $MarkerGroup{$id};
		$MarkerGroup{$id}=[$name,\@dat];
	}
}
close I;
ddx \%MarkerGroup;

open T,'<',$mapf or die $!;
<T>;
die "[x]Not a Tmap output file !" unless /^name/;
while (<T>) {
	my ($id,$cm)=split /\s+/;	# the file is TSV, but, ID comes with tailing 20h ...
	my ($name,$datref)=@{$MarkerGroup{$id}};
	push @MarkerOrder,$name;
}
