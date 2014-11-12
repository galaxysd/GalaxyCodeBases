#!/usr/bin/env perl
use strict;
use warnings;

my $minRefLen = 50;
my $maxRefLen = 100;
($minRefLen,$maxRefLen) = (40,220);

open I,'<','hg19chr18.bed.frag.trim' or die;

my $lines = 0;
while (<I>) {
	chomp;
	my ($chr,$s,$e,$len) = split /\t/;
	next if $len < $minRefLen or $len > $maxRefLen;
	++$lines;
}
my $each = int($lines/10);

print "Lines of [$minRefLen,$maxRefLen]: $lines\nEach of 10 win: $each\n\n";

seek I,0,0;

my ($len0,$s0) = 0;
while ($len0 < $minRefLen or $len0 > $maxRefLen) {
	$_=<I>;
	(undef,$s0,undef,$len0) = split /\t/;
}
my $i=1;

my @Rate = (
	[.1,.1],[1,0],[0,1],[.5,.1],[.1,.5],
	[.7,.3],[.3,.7],[.1,.5],[.5,.1],[.5,.5]
);
my $k=0;

open O,'>','zone.lst' or die;
my ($chr,$s,$e,$len);
while (<I>) {
	chomp;
	($chr,$s,$e,$len) = split /\t/;
	next if $len < $minRefLen or $len > $maxRefLen;
	++$i;
	unless ($i % $each) {
		print join("\t",$k,$s0,$e,@{$Rate[$k]}),"\n";
		print O join("\t",$k,$s0,$e,@{$Rate[$k]}),"\n";
		do {
			$_=<I>;
			last unless defined $_;
			(undef,$s0,undef,$len) = split /\t/;
		} while ($len < $minRefLen or $len > $maxRefLen);
		++$i;
		++$k;
	}
}
$k=9;
print join("\t",$k,$s0,$e,@{$Rate[$k]}),"\n";
print O join("\t",$k,$s0,$e,@{$Rate[$k]}),"\n";

close I;
close O;

__END__
scp sim0.pl gaoshengjie@192.168.9.100:/ifs4/BC_CANCER/PROJECT/qbe130308_gaosj_Methylation/sim
