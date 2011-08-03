#!/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 0.1.0 @ 20110803
=cut
use strict;
use warnings;

die "Usage: $0 <datfile> [max-depth-to-plot]\n" unless @ARGV;
my $name=$ARGV[0];
my $maxdep=$ARGV[1];

my $WinSize=0;
my (%GCdat,@GCvalues);
open IN,'<',$name or die "Error openimg $name: $!\n";
while(<IN>) {
    if (/^#WinSize=(\d+)/) {
        $WinSize=$1;
        next;
    }
    next if /^(#|$)/;
    chomp;
    my ($gc,undef,$depcnt,$mean)=split /\s+/;
    $gc -= 0.5;
    $GCdat{$gc}=[$depcnt,$mean];
    #print "{$gc}=[$depcnt,$mean]\n";
    push @GCvalues,$gc;
}
close IN;

open OUT,'>',$name.'.gc' or die "Error openimg ${name}.gc: $!\n";
print OUT "#WinSize: $WinSize\n",join("\t",'#GC%(left border)','MeanDepth'),"\n";
my $step=5; # 5%
my $lastGC=0;
my ($summer,$counter)=(0,0);
for my $gc (0..100) {
    if ($gc-$lastGC>=$step) {
        print OUT join("\t",$lastGC,$summer/$counter),"\n" if $counter;
        $lastGC=$gc;
        ($summer,$counter)=(0,0);
    }
    next unless exists $GCdat{$gc};
    $summer += $GCdat{$gc}->[1] * $GCdat{$gc}->[0];
    $counter += $GCdat{$gc}->[0];
}
#print OUT join("\t",$GCvalues[-1],$summer/$counter),"\n" if $counter;
close OUT;

chomp(my $user=`id -nru`);
my ($gnuplot,$font);
if ($user eq 'galaxy') {
    $gnuplot='gnuplot';
    $font='/home/galaxy/fonts/arial.ttf';
} else {
    $gnuplot='/opt/blc/genome/biosoft/gnuplot-4.4.0/bin/gnuplot';
    $font='/ifs1/ST_ASMB/USER/huxuesong/public/fonts/arial.ttf';
}
open P,'>',$name.'.dem' or die "Error openimg $name.dem: $!\n";
my $yrange='#set yrange [0:10]';
if ($maxdep) {
	$yrange="set yrange [0:$maxdep]";
}
print P <<"CODE";
reset
set xrange [-1:101]
$yrange
#set y2range [0:100]
set border 11
set xtics 0,5 nomirror
set ytics nomirror
set y2tics nomirror

set term png font "$font" 16 size 1600,1200
set output "$name.png"
set title "GC-dep plot of window_size $WinSize"
set xlabel "GC %"
set ylabel "Depth"
set y2label "Ref. Windows Count"
set boxwidth 0.8
set style fill empty
plot '$name' using 1:6:5:9:8 with candlesticks lt 3 lw 2 title 'depth' whiskerbars 0.5, \\
     ''      using 1:7:7:7:7 with candlesticks lt -1 lw 2 notitle, \\
     ''      using 1:4 lt rgb "red" lw 2 title 'mean depth', \\
     ''      using 1:3 with points axis x1y2 lt rgb "brown" lw 1 notitle, \\
     ''      using 1:11 with points lt rgb "royalblue" lw 1 title 'max depth', \\
     ''      using 1:3 smooth csplines axis x1y2 lt rgb "brown" lw 1 title 'Ref. win. cnt.'
#pause -1 "Hit return to continue"
CODE
close P;

system($gnuplot,$name.'.dem');
__END__

