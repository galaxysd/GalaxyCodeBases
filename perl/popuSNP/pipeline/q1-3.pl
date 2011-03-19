#!/bin/env perl
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.1;
my $SCRIPTS='/panfs/GAG/huxuesong/scripts';

our $opts='i:o:s:l:c:bv';
our($opt_i, $opt_s, $opt_l, $opt_o, $opt_v, $opt_b, $opt_c);

our $desc='1.filter fq, 2.soap, 3.rmdup';
our $help=<<EOH;
\t-i FASTAQ file path (./0rawfq) with *.adapter.list
\t-s Sample list (sample.lst) in format: /^Sample\\tLib\$/
\t-l Lib iNFOrmation (lib.nfo) in format:
\t   /^Lib\\tMax_Readlen\\tInsize_min\\tInsize_max\$/
\t-c Chromosome length list (chr.len) in format: /^ChrName\\s+ChrLen\\s?.*\$/
\t-o Project output path (.), will mkdir if not exist
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_i='./0rawfq' if ! $opt_i;
$opt_s='sample.lst' if ! $opt_s;
$opt_l='lib.nfo' if ! $opt_l;
$opt_o='.' if ! $opt_o;
$opt_c='chr.len' if ! $opt_c;

print STDERR "From [$opt_i] to [$opt_o] refer to [$opt_s][$opt_l][$opt_c]\n";
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
my (%SampleLib,%LibSample,%SampleMaxReadLen,%LibInsSize,%ChrLen);
open SAMPLE,'<',$opt_s or die "Error opening $opt_s: $!\n";
while (<SAMPLE>) {
	chomp;
	s/\r//g;
	my ($sample,$lib)=split /\t/;
	#print "Sample[$sample]\tLib[$lib]\n" if $opt_v;
	push @{$SampleLib{$sample}},$lib;
	#$LibSample{$lib}=$sample;
	push @{$LibSample{$lib}},$sample;
	$SampleMaxReadLen{$sample}=0;
}
close SAMPLE;
open LIBNFO,'<',$opt_l or die "Error opening $opt_l: $!\n";
while (<LIBNFO>) {
	chomp;
	s/\r//g;
	my ($lib,$maxrlen,$insmin,$insmax)=split /\s+/;
	my ($sample,$t);
	$t = ($sample) = @{$LibSample{$lib}};
	die "[!]Sample with Lib conflict($t) for Lib:[$lib]. Check [$opt_s] !\n" if $t != 1;	# Check if a lib -> more samples
	$SampleMaxReadLen{$sample} = $maxrlen if $SampleMaxReadLen{$sample} < $maxrlen;
	$LibInsSize{$lib}=[$insmin,$insmax];
}
close LIBNFO;
if ($opt_v) {
	for my $k (sort keys %SampleLib) {
		print "Sample:[$k]\tMaxRLen:[$SampleMaxReadLen{$k}]\tLib:[",join('],[',@{$SampleLib{$k}}),"]\n"
	}
	for my $k (sort keys %LibInsSize) {
		print "Lib:[$k]\tInsize:[",join(',',@{$LibInsSize{$k}}),"]\n";
	}
}
open CHRLEN,'<',$opt_c or die "Error opening $opt_c: $!\n";
while (<CHRLEN>) {
	chomp;
	s/\r//g;
	my ($chr,$len)=split /\s+/;
	#print "Chr:[$chr]\tLen:[$len]\n" if $opt_v;
	$ChrLen{$chr}=$len;
}
close CHRLEN;
if ($opt_v) {
	for my $k (sort keys %ChrLen) {
		print "Chr:[$k]\tLen:[$ChrLen{$k}]\n"
	}
}
=pod
# Check if a lib -> more samples
my %t;
++$t{$_} for (values %LibSample);
my $x=keys %t;
my $y=keys %SampleLib;
die '[!]',$y-$x," sample(s) with lib conflict. Check [$opt_s] !\n" if $x != $y;
=cut



system('mkdir','-p',$opt_o);


#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
#print "\nPlease use the following command to batch qsub:\033[32;1m
#find ${opt_o}/ -name '*.sh' | while read ll; do qsub -l vf=2G -cwd \$ll; done\n\033[0;0m\n";
__END__
