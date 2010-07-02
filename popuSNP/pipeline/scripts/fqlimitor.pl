#!/bin/env perl
use strict;
use warnings;

unless (@ARGV){
	print "perl $0 <BP cutoff> <in path/> <out path/> <.ext> <fq1 mainame> [fq2 mainame]\n";
	exit;
}

my ($maxBP,$inP,$outP,$ext,$fq1,$fq2) = @ARGV;
my $PE=1;
unless ($fq2) {
	$PE=0;
	#$maxBP /= 2;
}
my ($FQA,$FQB);

if ($ext eq '.fq') {
	open $FQA,'<',$inP.$fq1.$ext;
	open $FQB,'<',$inP.$fq2.$ext if $PE;
} elsif ($ext =~ /\.gz$/) {
	open $FQA,'-|',"gzip -dc ${inP}${fq1}${ext}";
	open $FQB,'-|',"gzip -dc ${inP}${fq2}${ext}" if $PE;
}
open OUTA,'>',$outP.$fq1.$ext;
open OUTB,'>',$outP.$fq2.$ext if $PE;

#my ($maxreadlenA,$copyedreadsA,$inbpA,$outbpA,$readlenA)=(0);
#my ($maxreadlenB,$copyedreadsB,$inbpB,$outbpB,$readlenB)=(0);
#my ($read_numA,$read_numB) = (0,0);

sub readfq($$) {
	my ($FH,$Aref)=@_;
	my ($readlen,$maxreadlen,$copyedreads,$bp,$Res)=@$Aref;
	my $line=<$FH>;
	return undef unless defined $line;
	$Res=$line;
	chomp($line=<$FH>);
	$readlen=length $line;
	$maxreadlen = $readlen if $maxreadlen < $readlen;
	$bp += $readlen;
	$Res .= $line."\n";
	$line=<$FH>; $Res .= $line;
	$line=<$FH>; $Res .= $line;
	++$copyedreads;
	return [$readlen,$maxreadlen,$copyedreads,$bp,$Res];
}

my ($DatRefA,$DatRefB)=([0,0,0,0,''],[0,0,0,0,'']);
while ($maxBP > 0) {
	$DatRefA=&readfq($FQA,$DatRefA);
	print OUTA pop @$DatRefA;
	$maxBP -= $$DatRefA[0];
	if ($PE) {
		$DatRefB=&readfq($FQB,$DatRefB);
		print OUTB pop @$DatRefB;
		$maxBP -= $$DatRefB[0];
	}
}
close $FQA;
close OUTA;
if ($PE) {
	close $FQB;
	close OUTB;
}

sub donfo($$) {
	my ($fq,$Aref)=@_;
	open NFO,'>',$outP.$fq.'.nfo';
	my $t=$inP.$fq.$ext;
# ($readlen,$maxreadlen,$copyedreads,$bp,$Res)=@$Aref;
	print NFO "$$Aref[2] parsed in [$t]
Using [NULL], No Adapter List used.
MaxReadLen\t$$Aref[1]
InReads\t$$Aref[2]
InBPs\t$$Aref[3]
FiltedReads\t0
CopyedReads\t$$Aref[2]
OutReads\t$$Aref[2]
OutBP\t$$Aref[3]
All done !\n";
close NFO;
}

&donfo($fq1,$DatRefA);
&donfo($fq2,$DatRefB) if $PE;
