#!/usr/bin/perl -w
#
use strict;

my $dataFile = $ARGV[0];
my $out_dir = $ARGV[1];
my $xrange = $ARGV[2];
my $xtic = $ARGV[3];
my $yrange = $ARGV[4];
my $ytic = $ARGV[5];
my $pro_name = $ARGV[6];
my $title = "Sequencing depth distribution";
my $font = 'font "Arial-Bold,24"';
my $titleFont = 'font "Arial-Bold,32"';
my $legendSize = 28;
my $xScale = 1;
my $yScale = 1;
my $pointSize = 1.3;



my $script = << "script";
set terminal postscript eps enhanced defaultplex leveldefault color colortext dashed dashlength 3.0 linewidth 1.0 butt palfuncparam 2000,0.003 "Helvetica" 16 

#set title '$title' $titleFont
#set lmargin 10
#set bmargin 7
set yrange [0:$yrange]
#set y2range [0:60]
set xrange [-0.5:$xrange]
set ytics $ytic $font
#set y2tics 10 $font
set xtics $xtic $font
set ytics nomirror

set boxwidth 0.8 absolute
#set style fill   solid 1.00 border  0
set style histogram clustered gap 1 title  offset character 0, 0, 0
set datafile missing '-'
set style data histograms
set xtics border in scale 1,0.5 nomirror offset character -1, -0.5, 10 
#set xtics border in scale 1,0.5 nomirror  offset character 0, -0.5, 10 

set style line 3  linetype 1 linecolor rgb "navy"   linewidth 5.000 pointtype 0 pointsize $pointSize
set style line 2  linetype 1 linecolor rgb "green"  linewidth 5.000 pointtype 0 pointsize $pointSize
set style line 1  linetype 1 linecolor rgb "pink"    linewidth 5.000 pointtype 11 pointsize $pointSize
set style line 4  linetype 1 linecolor rgb "red"    linewidth 5.000 pointtype 11 pointsize $pointSize

#set style histogram rowstacked


set key horizontal Left reverse enhanced autotitles columnhead nobox spacing 1.5 width 4.4 top right 
#unset key
#set key outside right spacing 1.2
set size $xScale,$yScale

set origin 0,0
set output "$out_dir/$pro_name.depthdistribution.eps"
set xlabe "\\nDepth(X)\\n\\n" $font
#set xlabe "\\n16mer(X)\\n\\n" $font
set ylabe "Percentage(\%)\\n" $font
#set y2labe "\\nGC(\%)\\n" $font


plot '$dataFile' using 1:3 title  "{/Arial=$legendSize $pro_name}" with linespoints ls 3 smooth csplines, '' u 1:4 title  "{/Arial=$legendSize Poisson}" with linespoints ls 2 smooth csplines


script
print "$script\n";
