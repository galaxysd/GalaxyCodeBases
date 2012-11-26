#!/bin/env perl
use strict;
use warnings;

die "Usage: $0 <in bam> <out path>\n" if @ARGV<2;

my $in=shift;
my $outp=shift;
my $infilename = (split /\//,$in)[-1];
$infilename =~ s/\.bam$//i;

my $SAMTOOLS='samtools';
system($SAMTOOLS) == -1 and $SAMTOOLS='/usr/local/bin/samtools';

my %Chr2ID;
open H, '-|', "$SAMTOOLS view -H $in" or die "Error opening $in: $!\n";
my @HEADS = <H>;
close H;
for (@HEADS) {
	my ($type,$chr,$len)=split /\t/;
	next unless $type eq '@SQ';
	$chr=(split(':',$chr))[1];
	my $str = $chr;
	if ($str =~ /gb\|([\w.]+)\|/) {
		$str=$1;
		$str=~s/\.\d+$//;
	} else {
		$str =~ s/\W/-/g;
		$str=~s/-+$//;
	}
	print "$str, $chr, $len";
	die if exists $Chr2ID{$chr};
	$Chr2ID{$chr} = $str;
}

my %Chr2FH;
for my $k (keys %Chr2ID) {
	my $out = "$outp/$infilename.$Chr2ID{$k}.bam";
	open $Chr2FH{$k},'|-', "$SAMTOOLS view -bS - >$out" or die "Error opening $out: $!\n";
	print { $Chr2FH{$k} } @HEADS;
}
open I, '-|', "$SAMTOOLS view $in" or die "Error opening $in: $!\n";
while (<I>) {
	my @items = split /\t/;
	next if $items[2] eq '*';
	my $fh = $Chr2FH{$items[2]} or die;
	print $fh $_;
}
close I;
close $Chr2FH{$_} for keys %Chr2FH;
__END__

