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
my $tt;
testundef($tt);

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

my ($FileCount,$MaxLen,%ChrLen,%IFH,%IFname,%CVGDAT,%CVGPOS)=(0,0);
open NFO,'<',$nfo or die "Error opening $nfo: $!\n";
print STDERR "\nLoading Chr Len ...\n[";
while (<NFO>) {
	my ($id,$len) = split /\s+/;
	$ChrLen{$id}=$len;
	$MaxLen = $len if $MaxLen < $len;
	print STDERR "$id:$len   ";
}
close NFO;
warn "\b\b\b].\n";
#ddx \%ChrLen;
open LST,'<',$in or die "Error opening $in: $!\n";
print STDERR "\nLoading Coverage files ...\n[";
while (<LST>) {
	my ($id) = split /\s+/;
	next if exists $IFname{$id};
	$IFname{$id}="./cvg/$id/cvg_$id.depth.xz";
	my $size = -s $IFname{$id};
#print "$IFname{$id}\n";
	next unless $size;
	$IFH{$id} = &openfile($IFname{$id});
	$CVGDAT{$id} = memalloc($MaxLen);
	#$CVGPOS{$id} = 0;
	print STDERR "$id:$size, ";
	++$FileCount;
}
close LST;
warn "\b\b].\n$FileCount file(s) opened.\n";

warn "\nReading Coverage data of $FileCount file(s) ...\n";
my $CurrentChrID='';
my ($fh,$chrid,$tmp);

for my $id (sort keys %IFH) {
	$fh = $IFH{$id};
	($tmp = <$fh>) =~ /^>(\S+)/;
	$chrid = $1;
	#while ( ($tmp = <$fh>) =~ /^>(\S+)/ ) { $chrid = $1; }; # How to push back ?
	unless (defined $tmp) {	# file end
		close $fh;
		delete $IFH{$id};
	}
	$CurrentChrID = $chrid if $CurrentChrID eq '';
	die "[x]ChrID order not match for [$chrid] in [$IFname{$id}] & [$CurrentChrID] !" if $CurrentChrID ne $chrid;
}

while (%IFH) {
	my $fhCount = (keys %IFH);
	die "[x]Depth file not match for $fhCount < $FileCount\n" if ($fhCount != $FileCount);
	print STDERR ">$CurrentChrID:     ";
	my $Filechrid = $CurrentChrID;
	$CurrentChrID = '';
	my $cnt=1;
	for my $id (sort keys %IFH) {
		$fh = $IFH{$id};
		$CVGPOS{$id}=0;
		printf STDERR "\b\b\b\b%4d",$cnt;
		while ( ($tmp = <$fh>) ) {
			last if $tmp =~ /^>/;
			my @dat = split /\s+/,$tmp;	# no need to chomp for /\s+/
			for my $v (@dat) {
				#pushdepvalue($CVGDAT{$Filechrid},$CVGPOS{$id},$v);
				++$CVGPOS{$id};
			}
#print STDERR "\n${Filechrid}[$id]$CVGPOS{$id}<$tmp>\n";
		}
		die "[x]Depth file format error in Chr:$chrid for $CVGPOS{$id} != $ChrLen{$Filechrid}." if $CVGPOS{$id} != $ChrLen{$Filechrid};
		unless (defined $tmp) {	# file end
			close $fh;
			delete $IFH{$id};
			next;
		}
		$tmp =~ /^>(\S+)/;
		$chrid = $1;
		$CurrentChrID = $chrid if $CurrentChrID eq '';
		die "\n[x]ChrID order not match for [$chrid] in [$IFname{$id}] & [$CurrentChrID] !" if $CurrentChrID ne $chrid;
		++$cnt;
	}
	warn "\b\b\b\bdone !\n";
}

for my $id (keys %CVGDAT) {
	memfree($CVGDAT{$id});
}
__END__



__C__

#include <stdlib.h>
#include <stdio.h>

void* memalloc(size_t size) {
	void * pt = malloc(size * 2);
	if (pt == NULL) {
		fputs("[x]Not enough memory !\n",stderr);
		exit(1); // Is it safe to exit here ?
	}
	return pt;
}
void* memrealloc(void *ptr, size_t size) {
	void * pt = realloc(ptr, size * 2);
	if (pt == NULL) {
		fputs("[x]Not enough memory !\n",stderr);
		free(ptr);
		exit(1); // Is it safe to exit here ?
	}
	return pt;
}
void memfree(void* pt) {
	free(pt);
}
void pushdepvalue(void* pt, size_t position, uint16_t value) {
	printf("@%zx:%d\t",position,value);
	*((uint16_t*)pt + position) = value;
}


void* testundef(void *ptr) {
	if (ptr == NULL) {
		fputs("It is NULL.\n",stderr);
	} else {
		printf("[%Lx]\n",ptr);
	}
}
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
