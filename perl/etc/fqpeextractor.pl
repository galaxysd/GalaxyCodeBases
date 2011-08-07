#!/bin/env perl
use strict;
use warnings;

die "Usage: $0 <len> <pairs> <fq1> <fq2>\n" if @ARGV != 4;
my ($length,$pairs,$fq1,$fq2)=@ARGV;
my ($out1,$out2);
$out1=$fq1.'.out';
$out2=$fq2.'.out';
warn "Len=$length Read_Pairs=$pairs\nFQ1=$fq1\nFQ2=$fq2\n";

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

sub readfq($) {
	my $fh=$_[0];
	chomp(my $a=<$fh>) or return [];
	chomp(my $b=<$fh>) or return [];
	chomp(my $c=<$fh>) or return [];
	chomp(my $d=<$fh>) or return [];
	return [$a,$b,$c,$d];
}

my $FQa=openfile($fq1);
my $FQb=openfile($fq2);

open Oa,'>',$out1 or die "Error: $!\n";
open Ob,'>',$out2 or die "Error: $!\n";

my $i=0;
while($i<$pairs) {
	my $refA=readfq($FQa);
	my $refB=readfq($FQb);
	next if $$refA[1]=~/N/i;
	next if $$refB[1]=~/N/i;
	$$refA[1]=substr $$refA[1],0,$length;
	$$refA[3]=substr $$refA[3],0,$length;
	$$refB[1]=substr $$refB[1],0,$length;
	$$refB[3]=substr $$refB[3],0,$length;
	print Oa join("\n",@$refA),"\n";
	print Ob join("\n",@$refB),"\n";
	++$i;
}
close Oa;
close Ob;
close $FQa;
close $FQb;
__END__
