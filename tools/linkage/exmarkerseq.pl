#!/bin/env perl
use strict;
use warnings;

unless (@ARGV > 0) {
    print "perl $0 <fabychr_path> <marker_len> <marker_file> <out_file>\n";
    exit 0;
}

my ($fabychr_path,$marker_len,$marker_file,$out_file)=@ARGV;
$marker_len=int $marker_len;
die "[x]marker_len must > 10 !\n" if $marker_len < 10;
open M,'<',$marker_file or die $!;
<M>;
my $t=tell M;
$_=<M>;
my $ChrID=(split /\t/)[0];
seek M,$t,0;

my $seq;
open F,'<',$fabychr_path."/$ChrID.fa" or die $!;
while(<F>) {
	chomp;
	my ($id)=split / /,$_,2;
	$id =~ s/^>//;
	die "[x]Wrong RefSeq ! ($ChrID != $id)\n" if $ChrID ne $id;
	$/=">";
	$seq=<F>;
	chomp $seq;
	$seq =~ s/\s//g;
	$/="\n";
}
close F;

open O,'>',$out_file or die $!;
while (<M>) {
	my $pos=(split /\t/)[1];
	my ($left,$right)=($pos-$marker_len,$pos+$marker_len);
	my $marker_seq=substr($seq,$left-1,$marker_len).'N'.substr($seq,$pos,$marker_len);
	print O ">${ChrID}_${pos}m\n$marker_seq\n";
}
close M;
close O;

__END__
./exmarkerseq.pl ./fabychr/ 32 ./dat20101214/Chr01.marker eChr01.marker
cat chrorder | while read a; do ./exmarkerseq.pl ./fabychr/ 32 ./dat20101214/${a}.marker ex32${a}.marker;done &
