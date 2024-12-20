#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Data::Dump qw(ddx);
use Galaxy::IO;
use Galaxy::SeqTools;

die "Usage: $0 <file1> <col1> <file2> <col2> <out.prefix>\n" if @ARGV<4;
my $file1=shift;
my $col1=shift;
my $file2=shift;
my $col2=shift;
my $outp=shift;
--$col1;
--$col2;

my (%dat,@remain);
open IA,'<',$file1 or die $!;
while (<IA>) {
	chomp;
	my @t = split /\s+/;
	$dat{$t[$col1]} = \@t;
}
close IA;

my %used;
open O,'>',"$outp.out" or die $!;
open IB,'<',$file2 or die $!;
while (<IB>) {
	chomp;
	my @t = split /\s+/;
	if (exists $dat{$t[$col2]}) {
		print O join("\t",@{$dat{$t[$col2]}},@t),"\n";
		++$used{$t[$col2]};
	} else {
		push @remain,\@t;
	}
}
close IB;
close O;
open OA,'>',"$outp.1" or die $!;
open OB,'>',"$outp.2" or die $!;
for (keys %dat) {
	print OA join("\t",@{$dat{$_}}),"\n" unless exists $used{$_};
}
for (@remain) {
	print OB join("\t",@{$_}),"\n";
}
close OA;
close OB;

__END__
awk '{print $2,$3,$4,$7,$8,$13,$14,$15}' felCat4.gtf |grep -v chrUn > felCat4.lst
perl -lane '$F[17]=~s/[";]//g;print join("\t",@F[0,3,4,6,17]) if $F[2] eq "transcript"' P_tigris.gene.gtf > P_tigris.lst
grep -v UNKNOWN P_tigris.lst > tigris.lst

./dojoin.pl felCat4.lst 6 tigris.lst 5 tigris2cat

awk '{print $2}' felCat4.lst|uniq|wc
awk '{print $2}' <(sort -k2 tigris2cat.out) |uniq|wc

sort -k2 tigris2cat.out > tigris2cat.sort

awk '{print $2,$3,$4,$7,$8,$13,$14,$15}' hg19.ucsc |grep -v chrUn|grep -v _random|grep -v _hap > hg19.lst
./dojoin.pl hg19.lst 6 tigris.lst 5 tigris2hum
sort -k2 tigris2hum.out > tigris2hum.sort

awk '{print $2,$3,$4,$7,$8,$13,$14,$15}' canfam3.ucsc |grep -v chrUn > canfam3.lst
./dojoin.pl canfam3.lst 6 tigris.lst 5 tigris2dog
sort -k2 tigris2dog.out > tigris2dog.sort

sort -nk10 ../rec.pas > rec.pas.nk10

perl -lane 'BEGIN {my ($cnt,%a)=(0);} ++$cnt;$l=log($F[9])/log(10);$l=int($l*5)/5;++$a{$l}; END { $t=0;print "-" x 75; for (sort {$a<=>$b} keys %a) { $t += $a{$_}; print "$_\t$a{$_}\t$t\t",int(100000*$t/$cnt)/1000; } }' rec.pas.nk10 > rec.pas.nk10.ratio

