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

open O,'>',"$outp.out" or die $!;
open IB,'<',$file2 or die $!;
while (<IB>) {
	chomp;
	my @t = split /\s+/;
	if (exists $dat{$t[$col2]}) {
		print O join("\t",@{$dat{$t[$col2]}},@t),"\n";
	} else {
		push @remain,\@t;
	}
}
close IB;
close O;
open OA,'>',"$outp.1" or die $!;
open OB,'>',"$outp.2" or die $!;
for (keys %dat) {
	print OA join("\t",@{$dat{$_}}),"\n";
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
