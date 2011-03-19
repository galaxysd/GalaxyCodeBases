#!/usr/bin/perl -w
#
use strict;

my $dataFile = $ARGV[0];
my $yrange = $ARGV[1] || 30;
my $ytics = $ARGV[2] || 5;
my $y2range = $ARGV[3] || 60;
my $y2tics = $ARGV[4] || 10;

my $title = "Sequencing depth distribution";
my $font = 'font "Arial-Bold,24"';
my $titleFont = 'font "Arial-Bold,32"';
my $legendSize = 28;
my $xScale = 1.6;
my $yScale = 1;
my $pointSize = 1.3;


my $script = << "script";
set terminal postscript eps enhanced defaultplex leveldefault color colortext dashed dashlength 3.0 linewidth 1.0 butt palfuncparam 2000,0.003 "Helvetica" 16 

#set title '$title' $titleFont
#set lmargin 10
#set bmargin 7
set yrange [0:$yrange]
set y2range [0:$y2range]
#set xrange [0:22]
set ytics $ytics $font
set y2tics $y2tics $font
set xtics 2 $font
set ytics nomirror

set boxwidth 0.8 absolute
#set style fill   solid 1.00 border  0
set style histogram clustered gap 1 title  offset character 0, 0, 0
set datafile missing '-'
set style data histograms
set xtics border in scale 1,0.5 nomirror   rotate by  -60 offset character -1, -0.5, 10 
#set xtics border in scale 1,0.5 nomirror  offset character 0, -0.5, 10 

set style line 1  linetype 1 linecolor rgb "navy"   linewidth 5.000 pointtype 5 pointsize $pointSize
set style line 2  linetype 1 linecolor rgb "green"  linewidth 5.000 pointtype 7 pointsize $pointSize
set style line 3  linetype 1 linecolor rgb "pink"    linewidth 5.000 pointtype 11 pointsize $pointSize

#set style histogram rowstacked


set key horizontal Left reverse enhanced autotitles columnhead nobox spacing 2.5 width 5.4 top left
#set key outside right spacing 1.2
set size $xScale,$yScale

set origin 0,0
set output "$dataFile.eps"
set xlabe "\\n\\nChromosome\\n\\n" $font
set ylabe "Sequencing depth(X)\\n" $font
set y2labe "\\nGC(\%)\\n" $font

plot '$dataFile' using 4:xtic(1) title "{/Arial=$legendSize Mode depth}" fs pattern 0 ls 1,'' u 3 title  "{/Arial=$legendSize Mean Depth}" with linespoints ls 2, '' using 5:xtic(1) title "{/Arial=$legendSize GC(\%)}" with linespoints axes x1y2 ls 3


script
print "$script\n";
