#!/bin/env perl
use lib '/share/raid010/resequencing/soft/lib';
use lib 'E:\BGI\toGit\perlib\etc';
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;
#use GalaxyXS::ChromByte 1.02;
#use DBI;
#
$main::VERSION=0.0.1;

our $opts='i:o:d:bv';
our ($opt_i, $opt_o, $opt_v, $opt_b, $opt_d);

our $help=<<EOH;
\t-i Group list (./veen.lst)
\t-o Output Stat (stat.txt)
\t-d Details dump to (details.lst)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_i='./veen.lst' if ! defined $opt_i;
$opt_o='stat.txt' if ! $opt_o;
$opt_d='details.lst' if ! $opt_d;

print STDERR "From [$opt_i] to [$opt_o] [$opt_d]\n";
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my ($bit,$i,%Log2,%files)=(1,0);
open( L,'<',$opt_i) or die "[x]Error: $!\n";
while (<L>) {
	use integer;
	chomp;
	my ($id,$file)=split /\t/;
	$bit *= 2;
	$Log2{$bit}=$i;
	$files{$id}=$file;
	++$i;
	warn "[!]$i: $id -> $file\n";
}

my @BitsA=sort {$a <=> $b} keys %Log2;
--$bit;

my %FH;
for my $id (keys %files) {
	open $FH{$id},'<',$files{$id} or die "[x]Error: $!\n";
}


close $FH{$_} for keys %files;



