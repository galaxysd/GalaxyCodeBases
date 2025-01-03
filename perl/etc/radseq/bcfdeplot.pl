#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);
use Galaxy::IO;
#use Galaxy::SeqTools;
#use Data::Dumper;

my $WinSize=10000;
$WinSize=20;
my $RegionBegin=151586958;
my $RegionEnd=152939134;
($RegionBegin,$RegionEnd) = (151586900,152939200);
my $Regions='gi|753572091|ref|NC_018727.2|:'."$RegionBegin-$RegionEnd";

die "Usage: $0 <mpileup bcf> <tfam file> <Win Size> <out>\n" if @ARGV<3;
my $bcfs=shift;
my $tfamfs=shift;
$WinSize=shift;
my $outfs=shift;

my (@tfamSamples,%tfamSamplePheno,%inFamily,@CaseS,@ControlS);
open F,'<',$tfamfs or die $!;
while (<F>) {
	next if /^(#|$)/;
	chomp;
	my ($family,$ind,$P,$M,$sex,$pho) = split /\t/;
	if ($pho == 1 or $pho == 2) {
		$tfamSamplePheno{$ind} = $pho;	# disease phenotype (1=unaff/ctl, 2=aff/case, 0=miss)
		if ($pho == 2) {
			push @CaseS,$ind;
		} else {
			push @ControlS,$ind;
		}
	} else { die; }	# $pho can only be 1 or 2
	push @tfamSamples,$ind;
	push @{$inFamily{$family}},$ind;
}
close F;

my $cmd = 'bcftools query -f \'%CHROM\t%POS\t%REF\t%DP[\t%SAMPLE=%DP]\n\' -r \'' . "$Regions' -s " . join(',',@CaseS,@ControlS);
warn "$cmd $bcfs |less -S\n";
my (%Count,%Filled);
my $fh = openpipe($cmd,$bcfs);
while (<$fh>) {
	chomp;
	my ($chr,$pos,$ref,$adp,@sampleDat) = split /\t/;
	#print "$chr,$pos,$ref:@sampleDat\n";
	my $refPos = $pos - $RegionBegin;
	my $zoneID = int($refPos/$WinSize);
	if ($ref eq 'N') {
		++$Filled{$zoneID}->[0];
		$Count{$zoneID}->[0] += $adp;
		next;
	}
	++$Filled{$zoneID}->[1];
	++$Filled{$zoneID}->[2];
	for (@sampleDat) {
		my ($id,$sdp) = split /=/,$_;
		my $StaType = $tfamSamplePheno{$id};
		$Count{$zoneID}->[$StaType] += $sdp;
	}
}
close $fh;

open O,'>',$outfs or die $!;
print O "# ",join("\t",qw/ZoneID ZonePos N Control Case/),"\n";
#ddx \%Filled;
#ddx \%Count;
my $maxZoneID = int(($RegionEnd-$RegionBegin)/$WinSize);
my $zoneRealBeginAdd = int($RegionBegin/$WinSize);
for my $zoneID (0 .. $maxZoneID) {
	my @ResLine=qw(? ? ?);
	for my $i (0 .. 2) {
		if (defined $Filled{$zoneID}->[$i]) {
			$ResLine[$i] = $Count{$zoneID}->[$i] / $Filled{$zoneID}->[$i] if defined $Count{$zoneID}->[$i];
		}
	}
	print O join("\t",$zoneID,($zoneID + $zoneRealBeginAdd)*$WinSize,@ResLine),"\n";
	#ddx $Count{$zoneID};
}
close O;

my $PlotCmd = <<"PLOTCMD";
set term png font "/opt/arial.ttf" 24 size 2400,900 truecolor linewidth 2
set datafile missing "?"
set autoscale	# scale axes automatically
unset log	# remove any log-scaling
unset label	# remove any previous labels
set xtic auto	# set xtics automatically
set ytic auto	# set ytics automatically
set title "Chr B2"
set xlabel "Window Position"
set ylabel "Coverage"

set yrange [0:500]

set format x "%.0f"
#set style fill transparent solid 0.5 noborder

set output "${outfs}z.png"
set style data linespoints
set xrange [152040000:152070000]
plot "$outfs" using 2:5 title 'Case' lc rgb 'red',\\
     '' using 2:4 title 'Control' lc rgb 'navy'

set output "$outfs.png"
set style data dots
set xrange [151586900:152939200]
plot "$outfs" using 2:5 title 'Case: Red' lc rgb 'red',\\
     '' using 2:4 title 'Control: Navy' lc rgb 'navy'

set style data lines
set output "${outfs}_Case.png"
plot "$outfs" using 2:5 title 'Case' lc rgb 'navy'
set output "${outfs}_Control.png"
plot "$outfs" using 2:4 title 'Control' lc rgb 'navy'
PLOTCMD

system "gnuplot << PLOTCMD\n$PlotCmd\nPLOTCMD\n";

__END__
#./bcfdeplot.pl mpileup_20150402HKT165931.bcf outA13.tfam outA13.depth10k
./bcfdeplot.pl mpileup_20150402HKT165931.bcf outA13.tfam 20 outA13.depth20
./bcfdeplot.pl mpileup_20150402HKT165931.bcf outA13.tfam 100000 outA13.depth100k


bcftools query -f '%CHROM,%POS\t%DP[\t%SAMPLE=%DP]\n' mpileup_20150402HKT165931.bcf|les

http://stackoverflow.com/questions/8714355/bash-turning-multi-line-string-into-single-comma-separated
grep -P '\t2$' outA13.tfam|awk -vORS=, '{print $2}'| sed 's/,$/\n/'
grep -P '\t2$' outA13.tfam|awk '{print $2}'| paste -d, -s

bcftools query -f '%CHROM,%POS\t%DP[\t%SAMPLE=%DP]\n' -s `grep -P '\t2$' outA13.tfam|awk '{print $2}'| paste -d, -s` mpileup_20150402HKT165931.bcf|les

gnuplot << PLOTCMD

set term png font "/opt/arial.ttf" 24 size 2400,900 truecolor linewidth 2
set datafile missing "?"
set autoscale	# scale axes automatically
unset log	# remove any log-scaling
unset label	# remove any previous labels
set xtic auto	# set xtics automatically
set ytic auto	# set ytics automatically
set title "Chr B2"
set xlabel "Window Position"
set ylabel "Coverage"

set yrange [0:500]

set format x "%.0f"
#set style fill transparent solid 0.5 noborder

set output "outA13z.png"
set style data linespoints
set xrange [152040000:152070000]
plot "outA13.depth10k" using 2:5 title 'Case' lc rgb 'red',\
     '' using 2:4 title 'Control' lc rgb 'navy'

set output "outA13.png"
set style data dots
set xrange [151586900:152939200]
plot "outA13.depth10k" using 2:5 title 'Case: Red' lc rgb 'red',\
     '' using 2:4 title 'Control: Navy' lc rgb 'navy'

set style data lines
set output "outA13_Case.png"
plot "outA13.depth10k" using 2:5 title 'Case' lc rgb 'navy'
set output "outA13_Control.png"
plot "outA13.depth10k" using 2:4 title 'Control' lc rgb 'navy'

PLOTCMD
