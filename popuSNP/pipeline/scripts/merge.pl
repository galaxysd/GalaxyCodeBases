#!/bin/env perl
use strict;
use warnings;

unless (@ARGV){
	print "perl $0 <merge list> <output prefix>\n";
	exit;
}
# $out.se, $out.log
my ($listf,$outf)=@ARGV;
open LST,'<',$listf or die "[x]Error opening $listf: $!\n";
my @files=<LST>;
chomp @files;
close LST;
if ($#files == 0) {
	link $files[0],"$outf.sp" or symlink $files[0],"$outf.sp" or exit 1;
} else {}
open LOG,'>',"$outf.log";
print LOG "done !\n";
close LOG;
