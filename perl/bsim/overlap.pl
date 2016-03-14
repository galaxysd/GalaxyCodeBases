#!/usr/bin/env perl
use strict;
use warnings;

die "Usage: $0 <SimLst> <AnalyseRes> [Shift=10]\n" if @ARGV <2;

my $total=shift;
my $result=shift;
my $shift=shift;
$shift = 10 unless defined $shift;
open TT,$total or die $!;
open RT,$result or die $!;
my %hash;

while(<TT>){
	chop;
	my @a=split;
	#print $a[3]."\n";
	next unless(/\w/);
	for( my $kk=$a[3]-$shift;$kk<$a[3]+$shift;$kk++){
		$hash{$kk}=1;
	}

}
close TT;
while(<RT>){
	chomp;
	my @a=split;
	#print $a[2]."\n";
	if(/RefCut/){
	if(exists($hash{$a[2]})){
		print $_."\n";
	}
	}


}
close RT;

__END__
./overlap.pl simed.lst simVir4_analyseAll.txt > simgot.lst

./overlap.pl simed.lst simVir4_analyseAll.txt 3|sort|uniq|wc -l
	142
./overlap.pl simed.lst simVir4_analyseAll.txt 4|sort|uniq|wc -l
	178
./overlap.pl simed.lst simVir4_analyseAll.txt 5|sort|uniq|wc -l
	186
./overlap.pl simed.lst simVir4_analyseAll.txt 10|sort|uniq|wc -l
	193
$ ./overlap.pl simed.lst simVir4_analyseAll.txt 20|sort|uniq|wc -l
	200

grep chr18 simVir4_analyseAll.txt|sort|uniq|wc -l
	247
grep RefCut simVir4_analyseAll.txt|sort|uniq|wc -l
	277

grep \> simout_*.Ref.fa | sed 's/^simout_m//'|sed 's/.Ref.fa:>/\t/'|sed 's/Ref_/Ref:/g'|sed 's/Vir_/Vir:/'|sed 's/R_/R:/'|sed 's/_/ /g'|cat -n >simed.lst

