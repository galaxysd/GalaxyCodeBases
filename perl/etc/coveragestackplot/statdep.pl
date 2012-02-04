#!/bin/env perl
use strict;
use warnings;
#use Inline 'Info', 'Force', 'Noclean';
use Inline 'C';

use Data::Dump 'ddx';
#use Inline 'C' => 'Config' => 'OPTIMIZE' => '-march=core2 -mtune=generic -O2 -pipe' ;
=pod
my $text = join '', <>;           # slurp input file
my $vp = vowel_scan($text);       # call our function
$vp = sprintf("%03.1f", $vp * 100);  # format for printing
print "The letters in @ARGV are $vp% vowels.\n";
=cut
my $a=testm();
print testp($a);
print testf($a),"\n";

die "Usage: $0 <chr.info> <list> [output]\n" if @ARGV <2;
my ($nfo,$in,$out)=@ARGV;
$out=$in.'.depstat' unless $out;
warn "From [$in] with [$nfo] to [$out]\n";
die "[x]Empty file: [$nfo] !\n" unless -s $nfo;
die "[x]Empty file: [$in] !\n" unless -s $in;

sub openfile($) {
    my ($filename)=@_;
    my $infile;
    if ($filename=~/.xz$/) {
	    open( $infile,"-|","xz -dc $filename") or die "Error opening $filename: $!\n";
    } elsif ($filename=~/.gz$/) {
     	open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
    } elsif ($filename=~/.bz2$/) {
     	open( $infile,"-|","bzip2 -dc $filename") or die "Error opening $filename: $!\n";
    } else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
    return $infile;
}

my ($FileCount,%ChrLen,%IFH,%IFname,%CCGDAT)=(0);
open NFO,'<',$nfo or die "Error opening $nfo: $!\n";
print STDERR "\nLoading Chr Len ...\n[";
while (<NFO>) {
	my ($id,$len) = split /\s+/;
	$ChrLen{$id}=$len;
	print STDERR "$id:$len   ";
}
close NFO;
warn "\b\b\b].\n";
#ddx \%ChrLen;
open LST,'<',$in or die "Error opening $in: $!\n";
print STDERR "\nLoading Coverage files ...\n[";
while (<LST>) {
	my ($id) = split /\s+/;
	$IFname{$id}="./cvg/$id/cvg_$id.depth.xz";
	my $size = -s $IFname{$id};
	next unless $size;
	$IFH{$id} = &openfile($IFname{$id});
	print STDERR "$id:$size, ";
	++$FileCount;
}
close LST;
warn "\b\b].\n$FileCount file(s) opened.\n";

__END__



__C__

#include <stdlib.h>

void* testm() {
	void * pt = calloc(1, 65537);
	sprintf(pt,"[Test Strings:%c%c%c.]",33,9,64);
	return pt;
}
SV* testp(void* pt) {
	return (newSVpvf("<%s>\n", (char*)pt));
}
int testf(void* pt) {
	free(pt);
	return 0;
}

/*
// Find percentage of vowels to letters
double vowel_scan(char* str) {
	int letters = 0;
	int vowels = 0;
	int i = 0;
	char c;
	char normalize = 'a' ^ 'A';
	// normalize forces lower case in ASCII; upper in EBCDIC
	char A = normalize | 'a';
	char E = normalize | 'e';
	char I = normalize | 'i';
	char O = normalize | 'o';
	char U = normalize | 'u';
	char Z = normalize | 'z';

	while(c = str[i++]) {
		c |= normalize;
		if (c >= A && c <= Z) {
			 letters++;
			 if (c == A || c == E || c == I || c == O || c == U)
				 vowels++;
		}
	}

	return letters ? ((double) vowels / letters) : 0.0;
}
*/
