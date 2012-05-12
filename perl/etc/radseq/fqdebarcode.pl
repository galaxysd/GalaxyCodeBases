#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120330
=cut
use strict;
use warnings;
#use Data::Dump qw(ddx);

die "Usage: $0 <bar code> <fq1> <outprefix1>\n" if @ARGV<2;
my $bar=shift;
my $fq=shift;
my $out=shift;

my $Eseq="CTGCAG";
$Eseq = "TGCAG";
#my $EcutAt=5;

sub openfile($) {
    my ($filename)=@_;
    my $infile;
    if ($filename=~/.bz2$/) {
	    open( $infile,"-|","bzip2 -dc $filename") or die "Error opening $filename: $!\n";
    } elsif ($filename=~/.gz$/) {
     	open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
    } else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
    return $infile;
}

my $fqfh=openfile($fq);
my $prefix = $bar.$Eseq;
my $barLen = length $bar;
open O,'|-',"gzip -9c >$out.fq.gz" or die "Error opening $out.fq.gz with gzip: $!\n";
open L,'>',$out.'.log' or die "Error opening $out.log: $!\n";
select(L);
$|=1;
print L "Barcode:[$bar]($barLen), Enzyme:[$Eseq].\n => Prefix:[$prefix]\nFrom [$fq] to [$out.fq.gz]\n";

my ($fq1,$fq2,$fq3,$fq4);
my ($tot,$have,$nothave)=(0,0,0);
while (defined ($fq1=<$fqfh>)) {
	$fq2=<$fqfh>;
	$fq3=<$fqfh>;
	$fq4=<$fqfh>;
	if ($fq2 =~ /^$prefix/) {
		$fq2 = substr $fq2,$barLen;
		$fq4 = substr $fq4,$barLen;
		++$have;
	} else {
		$fq2 = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n";
		$fq4 = "################################\n";
		++$nothave;
	}
	++$tot;
	print O join('',$fq1,$fq2,$fq3,$fq4);
}

close $fqfh;
print L "Total Reads: $tot\nWith Prefix: $have , ",$have/$tot,"\nWithout pfx: $nothave , ",$nothave/$tot,"\n";
close L;
close O;
