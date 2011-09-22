#!/bin/env perl
use strict;
use warnings;

die "Usage: $0 <read1.fq.gz> <out>\n" if @ARGV != 2;
my ($fq1,$out)=@ARGV;

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
#my $FQb=openfile($fq2);

my (%ReadsinTile,%Xrange,%Yrange);

my $fqitem=&readfq($FQa);
while (@$fqitem > 0) {
    my $id=$$fqitem[0];
    my ($Tile,$X,$Y)=(split /[:#]/,$id)[2,3,4];
    #my ($Surface,$Swath,$subTile)=(split //,$Tile)[0,1,3];
    #print "$Surface,$Swath,$subTile\t$Tile,$X,$Y\n";
    ++$ReadsinTile{$Tile};
    ++$Xrange{$X}; ++$Yrange{$Y};
    $fqitem=&readfq($FQa);
}
close $FQa;

open OUT,'>',$out or die "Error opening $out: $!\n";
for my $tile (sort keys %ReadsinTile) {
    print OUT "Tile\t$tile\t$ReadsinTile{$tile}\n";
}
print OUT "X\t$_\t$Xrange{$_}\n" for sort {$a<=>$b} keys %Xrange;
print OUT "Y\t$_\t$Yrange{$_}\n" for sort {$a<=>$b} keys %Yrange;
close OUT;



