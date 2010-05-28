#!/bin/env perl
use strict;
use warnings;

unless (@ARGV){
	print "perl $0 <matrix> <faByChr path> <ReadLen> <sorted soap> <out prefix> <Diploid/Mono>\n";
	exit;
}
# $out.soap, $out.single, $out.log
my ($matrix,$ref,$readlen,$sp,$out,$DM)=@ARGV;
$DM = ($DM=~/^m/i)?'-m ':'';
#warn $DM;
my $bin='/nas/RD_09C/resequencing/soft/bin/SoapSNP';
$sp =~ /\/([^_\/]+)_([^_]+)\.sp$/;
my $sample=$1;
my $ChrID=$2.'.fa';

my $sh="$bin -i $sp -d $ref/$ChrID -o $out.glf -F 1 ${DM}-I $matrix -L $readlen 2>$out.log";
# if you need monoploid calling mode, -m must be just in front of -I
open OUT,'>',"${out}_glf.sh.archive" or warn "[!]Error opening ${out}_glf.sh.archive: $!\n";
print OUT "#!/bin/sh\n$sh\n";
close OUT;
system($sh) or system("echo done ! > $out.tag");
