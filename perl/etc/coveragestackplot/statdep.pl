#!/bin/env perl
use strict;
use warnings;
use Inline 'Info', 'Force', 'Noclean';
use Inline 'C';

#use Data::Dump 'ddx';
#use Inline 'C' => 'Config' => 'OPTIMIZE' => '-march=core2 -mtune=generic -O2 -pipe' ;

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
	$CVGDAT{$id} = chrdat_init($MaxLen);
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
	chrdat_resize($CVGDAT{$_},$ChrLen{$Filechrid}) for keys %IFH;
	for my $id (sort keys %IFH) {
		$fh = $IFH{$id};
		$CVGPOS{$id}=0;
		printf STDERR "\b\b\b\b%4d",$cnt;
		while ( ($tmp = <$fh>) ) {
			last if $tmp =~ /^>/;
			my @dat;
			for ( split /\s+/,$tmp ) {	# no need to chomp for /\s+/
				push @dat,$_ if /^\d+$/;
			}
			#next unless @dat;
			$CVGPOS{$id} = chrdat_push($CVGDAT{$id},@dat);
=pod
			for my $v (@dat) {
				#next unless $v =~ /^\d+$/;
#print STDERR "$id,$CVGPOS{$id},$v\t";
				pushValue($CVGDAT{$id},$CVGPOS{$id},$v);
				++$CVGPOS{$id};
			}
=cut
#print STDERR "\n${Filechrid}[$id]$CVGPOS{$id}<$tmp>\n";
		}
#my @t=chrdat_readzone($CVGDAT{$id},1,$ChrLen{$Filechrid});
#print STDERR "Read:",scalar @t,"[",join(',',@t),"]\n";
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
	warn "\b\b\b\bdone of $ChrLen{$Filechrid} !\n";
}

for my $id (keys %CVGDAT) {
	chrdat_free($CVGDAT{$id});
}
__END__



__C__

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

typedef struct __ChrDat {
	size_t position, len;	// position is [0,len-1]
	uint16_t * dat;
} ChrDat_t;

void *chrdat_init(size_t size) {
	ChrDat_t * pt = malloc(sizeof(ChrDat_t));
	pt->len = size;
	pt->position = 0;
	pt->dat = malloc(size * 2);
	if (pt->dat == NULL) {
		fputs("[x]Not enough memory !\n",stderr);
		free(pt);
		exit(1); // Is it safe to exit here ?
	}
	return pt;
}
void chrdat_free(void * ptr) {
	ChrDat_t* pt = ptr;
	free(pt->dat);
	free(pt);
}
void chrdat_resize(void * ptr, size_t size) {
	ChrDat_t* pt = ptr;
	void * newpt = realloc(pt->dat, size * 2);
	if (newpt == NULL) {
		fputs("[x]Not enough memory !\n",stderr);
		chrdat_free(pt);
		exit(1); // Is it safe to exit here ?
	}
	pt->dat = newpt;
	pt->len = size;
	pt->position = 0;	// rewind here.
}
void chrdat_rewind(void * ptr) {
	ChrDat_t* pt = ptr;
	pt->position = 0;
}
SV* chrdat_push(void * ptr, ...) {
	ChrDat_t* pt = ptr;
	uint16_t * pDat = pt->dat + pt->position;
	Inline_Stack_Vars;
	int i;
//fputs("-[",stderr);
	for (i = 1; i < Inline_Stack_Items; i++) {	// ST(0) is ptr
		*pDat++ = (uint16_t) SvIV(Inline_Stack_Item(i));
		++pt->position;
//fprintf(stderr, "%zd:%d\t", pt->position, *(pDat-1) );
		if (pt->position > pt->len) {
			fputs("[x]Array too long !\n",stderr);
			exit(2); // Is it safe to exit here ?
		}
	}
	//Inline_Stack_Void;
	Inline_Stack_Done;
//fprintf(stderr, "\b]-\n");
	//pt->position = pDat - pt->dat;
	return newSViv(pt->position);
}
void chrdat_readzone(void * ptr, size_t start, size_t end) {
	ChrDat_t* pt = ptr;
	if (start < 1 || end > pt->len || start > pt->len || end < 1) {
		fputs("[x]Out Range !\n",stderr);
		return NULL;
	}
	uint16_t * pDat = pt->dat + start -1;
	Inline_Stack_Vars;

	Inline_Stack_Reset;
	while (pDat < pt->dat + end) {
		Inline_Stack_Push(sv_2mortal(newSViv(*pDat++)));
	}
	Inline_Stack_Done;
}
