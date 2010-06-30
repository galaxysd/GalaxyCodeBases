#!/bin/env perl
use lib '/share/raid010/resequencing/soft/lib';
use lib 'E:/BGI/toGit/perlib/etc';
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.1;

our $opts='i:o:l:s:n:v:bd';
our($opt_i, $opt_o, $opt_v, $opt_b, $opt_d, $opt_l, $opt_s, $opt_n);

our $desc='';
our $help=<<EOH;
\t-i Genome FASTA file (./genome.fa)
\t-l Reads Length (90)
\t-n iNsert size (500)
\t-s Slip distance (1)
\t-o Output FQ prefix (./out)-N-RL-Sd_{1,2}.fa
\t-v Verbose level (undef=0)
\t-b No pause for batch runs
\t-d Debug Mode, for test only
EOH

ShowHelp();

$opt_i='./genome.fa' if ! $opt_i;
no warnings;
$opt_v=int $opt_v;
$opt_n=int $opt_n;
$opt_l=int $opt_l;
$opt_s=int $opt_s;
use warnings;
$opt_n=500 if ! $opt_n;
$opt_l=90 if ! $opt_l;
$opt_s=1 if ! $opt_s;
$opt_o='./out' if ! $opt_o;
my $outfile="${opt_o}-${opt_n}-${opt_l}-${opt_s}";
my $lenB=$opt_n-$opt_l;

my $maxN=int(.3*$opt_l);

print STDERR "From [$opt_i] to [${opt_o}]-${opt_n}-${opt_l}-${opt_s}_{1,2}.fa, [$opt_n][$opt_l][$opt_s]\n";
print STDERR "DEBUG Mode on !\n" if $opt_d;
print STDERR "Verbose Mode [$opt_v] !\n" if $opt_v;
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
open G,'<',$opt_i or die "Error opening $opt_i: $!\n";
open OA,'>',$outfile.'_1.fa' or die "Error opening ${outfile}_1.fa: $!\n";
open OB,'>',$outfile.'_2.fa' or die "Error opening ${outfile}_2.fa: $!\n";
my ($str,$a,$b,$t);
while (<G>) {
	s/^>//;
	my $title = $_;
	my $seqname = $1 if($title =~ /^(\S+)/);
	print STDERR "[!]Seq. >$seqname:\b";
	$/=">";
	my $seq=<G>;
	chomp $seq;
	$seq=~s/\s//g;
	$/="\n";
	my $Len=length $seq;
	#$seq='';
	print STDERR "\t$Len .\b";
	for (my $i=0;$i<$Len-$opt_n;$i+=$opt_s) {
		$str=substr $seq,$i,$opt_n;
		$a=substr $str,0,$opt_l;
		$t=$a=~tr/Nn/nN/;
		next if $t > $maxN;
		$b=substr $str,$lenB;
		$t=$b=~tr/Nn/nN/;
		next if $t > $maxN;
		print OA ">${seqname}_${i}_1 ",$i+1,"\n$a\n";
		print OB ">${seqname}_${i}_2 ",$i+$lenB+1,"\n$b\n";
	}
	warn "-\n";
}
close OA;
close OB;
close G;

#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__
