#!/bin/env perl
use strict;
use warnings;

die "Usage: $0 <ori-fq> <B-masked-fq> <outprefix>.(fq.xz|log)\n" if @ARGV <2;
my ($inA,$inB,$out)=@ARGV;
$out=$inA.'.merge' unless $out;
warn "From [$inA]&[$inB] to [$out].(fq.xz|log)\n";

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
sub getFQitem($) {
	my $fh=$_[0];
	my ($a,$b,$c,$d);
	defined($a=<$fh>) or return 0;
	defined($b=<$fh>) or return 0;
	defined($c=<$fh>) or return 0;
	defined($d=<$fh>) or return 0;
	return [$a,$b,$c,$d]
}
sub CnM($$) {
	my ($a,$b)=@_;
	my ($id,$seq)=($$a[0],$$a[1]);
	if ($id ne $$b[0]) {
		chomp $id;
		$id .= "\t".$$b[0];
	}
	if ($seq ne $$b[1]) {
		chomp $seq;
		$seq .= "\t".$$b[1];
	}
	return "${id}${seq}$$a[3]$$b[3]";
}

my ($Count1,$Count2,$CountPairs)=(0,0,0);
my $fhA=openfile($inA);
my $fhB=openfile($inB);
$|=1;
open(OUT,"|-","xz -c9e > \"${out}.fq.xz\"") or die "Error opening $out:$!\n";
my ($dat1,$dat2);
while (1) {
	$dat1=getFQitem($fhA);
	$dat2=getFQitem($fhB);
	if ($dat1 && $dat2) {
		print OUT CnM($dat1,$dat2);
		++$CountPairs;
	} elsif ($dat1) {
		print OUT join('',@$dat1);
		++$Count1;
	} elsif ($dat2) {
		print OUT join('',@$dat2);
		++$Count2;
	} else {
		last;
	}
}

close $fhA;
close $fhB;
close OUT;
open LOG,'>',"${out}.log" or die "Error opening $out.log:$!\n";
my $str = "Out Pairs: $CountPairs\nFQ1 over hang: $Count1\nFQ2 over hang:$Count2\n";
print $str;
print LOG "From [$inA]&[$inB] to [$out.fq.xz]\n$str";
close LOG;
#print $Count{1}+$Count{-1}," ,RF:$CountRF\t+$Count{1},-$Count{-1},z$Count{0}\t$Reads\n";
__END__
find . -name '*.r'|xargs -n1 ./samrstat.pl
find . -name '*.ro'|xargs -n1 cat|sort -n
