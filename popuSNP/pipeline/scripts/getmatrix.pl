#!/bin/env perl
use strict;
use warnings;

unless (@ARGV){
	print "perl $0 <megred.lst> <fabychr> <out prefix>\n";
	exit;
}
# $out.soap, $out.single, $out.log
my ($lst,$ref,$out)=@ARGV;
my $bin='/nas/RD_09C/resequencing/soft/bin/SoapSNP';

my ($maxRL,$maxSize,$maxFile,$maxChr)=(0,0);
open L,'<',$lst;
while (<L>) {
	chomp;
	my ($sample,$chr,$readlen,$file)=split /\t/;
#print "$sample,$chr,$readlen,$file\n";
	my $size = (-s $file) || 0;
	next if $readlen < $maxRL;	# the longer reads can stands for the shorter ones
	if ($size > $maxSize or $readlen > $maxRL) {
		$maxSize = $size;
		$maxFile = $file;
		$maxRL = $readlen;
		$maxChr = $chr;
	}
}
close L;
$maxRL = 101 if $maxRL == 0;
die "[x]No soap file found !\n" unless $maxFile;
print "$maxRL,$maxSize,$maxFile\n";

my $sh="$bin -i $maxFile -d $ref/$maxChr.fa -o /dev/null -F 1 -M $out.matrix0 -L $maxRL 2>$out.log";
open OUT,'>',"${out}_matrix.sh.archive" or warn "[!]Error opening ${out}_matrix.sh.archive: $!\n";
print OUT "#!/bin/sh\n$sh\n";
close OUT;
if ( system($sh) ) {
	warn "[!]system [$sh] failed: $?\n" and exit 1;
} else {
	system('mv','-f',"$out.matrix0","$out.matrix") and exit 0;	# the exit value is 0 for normal
	#exit 2;
}
