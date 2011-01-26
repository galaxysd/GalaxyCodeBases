#!/bin/env perl
use strict;
use warnings;
use lib '/nas/RD_09C/resequencing/soft/lib';
#use GalaxyXS::ChromByte;
use Data::Dump qw(ddx);

unless (@ARGV > 0) {
    print "perl $0 <group_file> <tmap output> <Out_file>\n";
    exit 0;
}

my ($in,$mapf,$outf)=@ARGV;
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
#ddx \%MarkerGroup;

open T,'<',$mapf or die $!;
$_=<T>;
die "[x]Not a Tmap output file !" unless /^name/;
open O,'>',$outf or die $!;
while (<T>) {
	my ($id,$cm)=split /\s+/;	# the file is TSV, but, ID comes with tailing 20h ...
	my ($name,$datref)=@{$MarkerGroup{$id}};
	#push @MarkerOrder,$id;
	#$MarkerCM{$id}=$cm;
	print O "#$name\n";
	for (@{$datref}) {
		print O $_,"\t$cm\n";
	}
	print O "#\n";
}

__END__
./groupreader.pl deChr12.group deChr12.imo deChr12.cm
cat ../chrorder | while read a;do ./groupreader.pl de$a.group de$a.imo de$a.cm;done &

http://www.maizemap.org/iMapDB/Overview/Anchoring_Rules.html
Anchoring Rules
Rules developed to help make unambiguous contig: genetic map assignments

Rule 1.
One contig should have one position on the integrated map.

Rule 2.
Contig:marker associations made by hybridization of a probe to a single BAC in the contig will be filtered out before anchoring. Contig:marker associations made by detection of a single BAC by PCR-based analysis of the BAC DNA pools will be accepted if no other conflicts exist for that contig.

Rule 3.
If two linked markers* hit one shared contig, the contig will be anchored to the position of those markers, even if the markers also hit other contigs.

Rule 4.
Multiple contigs can be assigned to a single position if they are detected by probes corresponding to closely linked markers.

Rule 5.
If a marker detects multiple BACs that are uniquely assembled in one contig, and no conflicts exist for that contig, the contig can be assigned to the locus corresponding to that marker.


*. The closely linked markers are those that are within 5 marker's positions or within 10 cM with their neighbors on the genetic map meanwhile the markers with same coordinates are considered as same position.
