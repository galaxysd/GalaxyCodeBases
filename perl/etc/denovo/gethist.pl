#!/bin/env perl
use strict;
use warnings;

die "Usage: $0 <col> <input[gzipped]> <output prefix>\n" if @ARGV < 3;
my ($col,$inf,$outf)=@ARGV;
warn "From Col=[$col] of [$inf] to [$outf]\n";

--$col;
my ($Sum,$Count,%Cnt)=(0,0);

sub openFH($) {
	my $inf = $_[0];
	my $FH;
	if ($inf =~ /\.gz$/i) {
		open $FH,'-|',"gzip -dc $inf" or die "Error opening $inf: $!\n";
	} elsif ($inf =~ /\.xz$/i) {
		open $FH,'-|',"xz -dc $inf" or die "Error opening $inf: $!\n";
	} elsif ($inf =~ /\.bz2$/i) {
		open $FH,'-|',"bzip2 -dc $inf" or die "Error opening $inf: $!\n";
	} else {
		open $FH,'<',$inf or die "Error opening $inf: $!\n";
	}
	return $FH;
}

my $FH = openFH($inf);
while (<$FH>) {
	my @Dat = split /\s+/;
	next if /^[#;]/;
	++$Cnt{$Dat[$col]};
	++$Count;
	$Sum += $Dat[$col];
}
close $FH;

open O,'>',"$outf.dat" or die $!;
print O "# TotalCount=$Count, ValueSum=$Sum, Average=",$Sum/$Count,"\n#Value\tCount\tCountRatio\tCumCountRatio\n";

my $Cum=0;
for my $k (sort { $a<=> $b } keys %Cnt) { # $Cnt{$a} <=> $Cnt{$b} || $a<=> $b
	$Cum += $Cnt{$k}/$Count;
	print O join("\t",$k,$Cnt{$k},$Cnt{$k}/$Count,$Cum),"\n";
}
close O;

__END__
set terminal png notransparent nointerlace size 1200,960 font '/opt/arial.ttf' 24
set output "plot.png"
set xlabel 'Depth'
set ylabel 'Count'
set xrange[2:462]
set yrange[0:320000]
plot 'tw.dat' using 1:2 with lines
