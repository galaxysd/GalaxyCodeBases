#!/bin/env perl
use strict;
#use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.1;

our $opts='i:o:lbv';
our($opt_i, $opt_o, $opt_l, $opt_v, $opt_b);

our $desc='Dynamic String Alignmenter';
our $help=<<EOH;
\t-i Input file for 2 strings divided by \\n (in.txt)
\t-l use local alignment instead of global one
\t-o output file (out.txt)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_i='in.txt' if ! $opt_i;
$opt_o='out.txt' if ! $opt_o;

print STDERR "From [$opt_i] to [$opt_o], \033[32;1m";
if ($opt_l) {
	print "LOCAL\033[0;0m alignment.\n";
} else {print "GLOBAL\033[0;0m alignment.\n"}
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
open IN,'<',$opt_i or die "[x]Error opening $opt_i: $!\n";
$a=<IN> or die "[x]Error: $opt_i is empty!\n";
$b=<IN> or die "[x]Error: $opt_i contains only 1 line!\n";
close IN;
chomp($a,$b);
if ($opt_v) {
	print "\033[32;1mStrings loaded:\nX: [\033[0;0m$a\033[32;1m]\nY: [\033[0;0m$b\033[32;1m]\033[0;0m\n";
}

my ($MATCH,$ChrMATCH)=(5,3);
my $INDEL=-6;	# INDEL Well, just linear gap score.
my %SCOR=(	# MISMATCH
	d => { d=>-1,o=>-3,s=>-2, },
	o => { o=>-2,s=>-3, },
	s => { s=>-2, },	# -2 for test, should be -1
);
sub scoring($$) {
	my ($a,$b)=@_;
	return $MATCH if $a eq $b;	# MATCH
	return $ChrMATCH if lc $a eq lc $b;	# Chr MATCH
	my $m = ($a =~ /\d/)?'d':(($a =~ /\w/)?'s':'o');
	my $n = ($b =~ /\d/)?'d':(($b =~ /\w/)?'s':'o');
	($m,$n)=sort ($m,$n);	# triangle is better than square.
#print "($m,$n)\n";
	return $SCOR{$m}{$n};	# MISMATCH
}

my @a=('',split //,$a);
my @b=('',split //,$b);
my (@matrix,@path);	# 2D array in perl can be malloced auto.ly
my ($i,$j,$sc);
my %COLOR=('1'=>'32','2'=>'33','3'=>'36','0'=>'0');
# Initialization
$matrix[0][0]=0;
$path[0][0]=0;	# 0 is better for local alignment to always stop @ 0=FALSE (?!)
for ($i=1;$i<=$#a;$i++) {$matrix[$i][0] = $i*$INDEL;$path[$i][0]=2;}
for ($j=1;$j<=$#b;$j++) {$matrix[0][$j] = $j*$INDEL;$path[0][$j]=3;}
# dynamic recursion
for ($i=1;$i<=$#a;$i++) {
	for ($j=1;$j<=$#b;$j++) {
		$matrix[$i][$j] = scoring($a[$i],$b[$j]) + $matrix[$i-1][$j-1];
		$path[$i][$j] = 1;	# case 0: i,j are aligned
		$sc=$matrix[$i-1][$j] + $INDEL;	# case 1: i aligned to -
		if ($sc > $matrix[$i][$j]) {
			$matrix[$i][$j] = $sc;
			$path[$i][$j] = 2;
		}
		$sc=$matrix[$i][$j-1] + $INDEL;	# case 2: j aligned to -
		if ($sc > $matrix[$i][$j]) {
			$matrix[$i][$j] = $sc;
			$path[$i][$j] = 3;
		}
		if ($opt_l and $matrix[$i][$j] < 0) {
			$matrix[$i][$j] = 0;
		}
	}
}
# Traceback (with shadow matrix)
my ($count,@ax,@ay,@as)=(0);	# first is 0.
--$i;--$j;	# outside for, $i will bigger than $#a by 1 step

while ($i>=0 and $j>=0) {
	if ($opt_l) {
	} else {
		$as[$count] = $matrix[$i][$j] unless $as[$count];
		#&traceback($a[$i],$b[$j],$path[$i][$j],$ax[$count],$ay[$count]);
		my $path=$path[$i][$j];
print "$i,$j,$count,$path,$matrix[$i][$j]\t";
		if ($path == 1) {
			push @{$ax[$count]},$a[$i];push @{$ay[$count]},$b[$j];
print "$a[$i],$b[$j] ";
			--$i;--$j;
		} elsif ($path == 2) {
			push @{$ax[$count]},$a[$i];push @{$ay[$count]},"-";
print "$a[$i],- ";
			--$i;
		} elsif ($path == 3) {
			push @{$ax[$count]},"-";push @{$ay[$count]},$b[$j];
print "-,$b[$j] ";
			--$j;
		} else {	# $path == 0
			@{$ax[$count]} = reverse @{$ax[$count]};
			@{$ay[$count]} = reverse @{$ay[$count]};
			++$count;--$i;--$j;
print "<$i,$j,$count,$path,$matrix[$i][$j]>\t";
		}
print "\n";
	}

}

if ($opt_v) {
	print "\033[32;1mScoring system:\033[0;0m\n$MATCH for match; $ChrMATCH for char in different case; $INDEL for gap\n";
	print "Mismatch between \\d as \033[32;1md\033[0;0m, [\\w^\\d] as \033[32;1ms\033[0;0m and others as \033[32;1mo\033[0;0m are:\n";
	for my $key (sort keys %SCOR) {
		for (sort keys %{$SCOR{$key}}) {
			print " ${key}-$_:\033[32;1m$SCOR{$key}{$_}\033[0;0m";
		}
		print "\n";
	}
	print "\033[32;1mDynamic matrix:\033[0;0m\n";
	print "         ";printf "  [%s]",$_ for (@b[1..$#b]);print "\n";
	for ($i=0;$i<=$#a;$i++) {
		if ($i==0) {print "    "}
		 else {printf " [%s]",$a[$i]}
		for ($j=0;$j<=$#b;$j++) {
			print "\033[",$COLOR{$path[$i][$j]},";1m";
			printf("%4d ",$matrix[$i][$j]);
			#print "\033[0;0m";
		}
		print "\033[0;0m\n";
	}
	print "\033[32;1mAlignment:\033[0;0m\n";
	for ($i=1;$i<=$count;$i++) {
		print "No. $count:\n";
		print 'X:[',@{$ax[$i-1]},"]\n",'Y:[',@{$ay[$i-1]},"]\n",'Score:',$as[$i-1],"\n";
	}
}

#END
my $stop_time = [gettimeofday];
print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__
