#!/bin/env perl
use strict;
use warnings;
use IO::Unread qw(unread);
use Data::Dump qw(ddx);

die "Usage: $0 <input>\n" if @ARGV < 1;
my ($inf)=@ARGV;

sub openfile($) {
	my ($filename)=@_;
	my $infile;
	if ($filename=~/.xz$/) {
			open( $infile,"-|","xz -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.gz$/) {
		open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.lz4$/) {
		open( $infile,"-|","/opt/bin/lz4c -dy $filename /dev/stdout") or die "Error opening $filename: $!\n";
	} else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
	return $infile;
}

my $IN = openfile($inf);
my $tmp = <$IN>;
unread $IN, $tmp;
my $str = (split /\t/,$tmp)[0];
my $len = length $str;
my $allcnt = 4 ** $len;

my ($Sum,$kMax,$kMin,$cMax,$cMin,%Hist)=(0,0,999999999,0,999999999);
my ($KmerSum,$K90,$K75,$K50,$K25,$K10,$Kadding)=(0,0,0,0,0,0,0);

while (<$IN>) {
	chomp;
	my ($str,$cnt) = split /\t/;
	#print "$str\t[$cnt]\n";
	++$Hist{$cnt};
	++$Sum;
	$KmerSum += $cnt;
	$kMax = $cnt if $cnt > $kMax;
	$kMin = $cnt if $cnt < $kMin;
}
close $IN;

for (values %Hist) {
	$cMax = $_ if $cMax < $_;
	$cMin = $_ if $cMin > $_;
}

open OUT,'>',"$inf.hist" or die;
print OUT "# K\t$len\t$allcnt
# KmerFreq_Min-Max\t$kMin\t$kMax
# Cnt_Min-Max\t$cMin\t$cMax
# CntOfFoundKmer\t$Sum\t",$Sum/$allcnt,"\n";

my @order = sort {$a <=> $b} keys %Hist;
for (@order) {
	$Kadding += $_ * $Hist{$_};
	if ( 10 * $Kadding <  $KmerSum ) {
		$K90 = $_;
	} elsif ( 4 * $Kadding <  $KmerSum ) {
		$K75 = $_;
	} elsif ( 2 * $Kadding <  $KmerSum ) {
		$K50 = $_;
	} elsif ( 4 * $Kadding <  $KmerSum * 3 ) {
		$K25 = $_;
	} elsif ( 10 * $Kadding <  $KmerSum * 9 ) {
		$K10 = $_;
	}
}
print OUT "# KmerSum: $KmerSum\n# K90: $K90, K75: $K75, K50: $K50, K25: $K25, K10: $K10\n";
my @t;
push @t, $_/$KmerSum for ($K90,$K75,$K50,$K25,$K10);
print OUT '# Kratio: ',join(', ',@t),"\n";

print OUT "# KmerFreq\tCntOfThisFreq\tKmerSumRatio\tCumulativeKmerSumRatio\n";
my $cumulative = 0;
for (@order) {
	$cumulative += $_*$Hist{$_};
	print OUT join("\t",$_,$Hist{$_},$_*$Hist{$_}/$KmerSum,$cumulative/$KmerSum),"\n";
}
close OUT;

__END__
cat kmerfreq.lst|while read a;do perl getkmhist.pl $a; done &

cat kmerfreq.lst|while read a;do echo perl "getkmhist.pl $a &"; done > t.sh
