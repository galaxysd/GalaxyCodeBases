#!/bin/env perl
use strict;
use warnings;

unless (@ARGV){
	print "perl $0 <sorted soap> <faByChr path> <ReadLen path> <out prefix>\n";
	exit;
}
# $out.soap, $out.single, $out.log
my ($sp,$ref,$readlen,$out)=@ARGV;
my $bin='/panfs/GAG/huxuesong/scripts/soapsnp';
$sp =~ /\/([^_\/]+)_([^_]+)\.sp$/;
my $sample=$1;
my $ChrID=$2.'.fa';
my @len=`find $readlen/$sample/ -name '*.ReadLen'`;
chomp @len;
$readlen=0;
for (@len) {
	open LEN,'<',$_ or die "[x]Error opening $_: $!\n";
	my $len = <LEN>;
	close LEN;
	chomp $len;
	$readlen=$len if $readlen<$len;
}
my $sh="$bin -i $sp -d $ref/$ChrID -o /dev/null -F 1 -M $out.matrix0 -L $readlen 2>$out.log";
open OUT,'>',"${out}_matrix.sh.archive" or warn "[!]Error opening ${out}_matrix.sh.archive: $!\n";
print OUT "#!/bin/sh\n$sh\n";
close OUT;
if ( system($sh) ) {
	system('mv','-f',"$out.matrix0","$out.matrix") or exit 1;
} else {
	system('mv','-f',"$out.matrix0","$out.matrix") or exit 1;	# the exit value is 0 for normal, I was planing to modify soapSNP ...
	warn "[!]system [$sh] failed: $?";
	exit 2;
}
