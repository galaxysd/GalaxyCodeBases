#!/usr/bin/perl -w
use strict;

if (@ARGV < 2) {
	die "Usage: $0 <divisor> <in> [out]\n";
}

$ARGV[2] = '-' if (@ARGV == 2);

open IN, $ARGV[1] or die $!;
open OUT, ">$ARGV[2]" or die $!;

while(<IN>){
	chomp;
	my @sp=split /\,/, $_;
	my $len=@sp;
	my @clum=(0,0,0,0,0,0,0,0,0,0,0,0);
	my @one=split /\=/,$sp[0];
	$one[0]=~/.*\"(.*)\".*/;
	$clum[0]=$1;
	$one[1]=~/\D+(\d+)\D+/;
#	print $1;exit;
	my $num=$1+1;
	$one[2]=~/\D+(\d+)/;
#	print $1;exit;
	my $last=$1/$ARGV[0];
	$clum[$num]=$last;
	my $n=1;
	while($n<$len){
		my @two=split /\=/,$sp[$n];
		$two[0]=~/\D+(\d+)\D+/;
		my $num=$1+1;
		#print $sp[$n];exit;
		#print $two[1];exit;
		$two[1]=~/\D+(\d+)/;
		#print $1;exit;
		my $last=$1/$ARGV[0];
		$clum[$num]=$last;
		$n++;
	}
	print OUT "$clum[0]\t$clum[1]\t$clum[2]\t$clum[3]\t$clum[4]\t$clum[5]\t$clum[6]\t$clum[7]\t$clum[8]\t$clum[9]\t$clum[10]\n";
}

close IN;
close OUT;

__END__
./bsuit check sj/j90.ini |tail -n11|head -n10 > x90.txt
./getchecked.pl 1 x90.txt
