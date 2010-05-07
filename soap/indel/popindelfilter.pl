#!/bin/env perl
use lib '/share/raid010/resequencing/resequencing/tmp/bin/annotation/glfsqlite';
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.2;

our $opts='i:o:s:m:x:q:l:bv';
our ($opt_i, $opt_o, $opt_s, $opt_v, $opt_b, $opt_m, $opt_x, $opt_q, $opt_l);

our $help=<<EOH;
\t-i PSNP list (./psnp.lst) for chrid.individual.finalSNPs
\t-l Indel list (./indel.lst) in [SampleID\\tpath/to/indel-result.list]
\t-s GLF list (./glf.lst), will use \$1 of (([^/]+)/[^/]+$) for sample names
\t-o Output Prefix (./indel_f/f_)
\t-m minimal depth of the indel (3)
\t-x maximal depth of the indel (20)
\t-q the minimal quality of the indel (20)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_o='./indel_f/f_' if ! defined $opt_o;
$opt_i='./psnp.lst' if ! $opt_i;
$opt_s='./glf.lst' if ! $opt_s;
$opt_m=3 if ! $opt_m;
$opt_x=20 if ! $opt_x;
$opt_q=20 if ! $opt_q;

$opt_i =~ s/\/$//;
print STDERR "From [$opt_i] to [$opt_o], with [$opt_s]($opt_m,$opt_x)[$opt_q]/\n";
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];

system('mkdir','-p',$opt_o);
system('rmdir',$opt_o) if $opt_o =~ m#/[\s.]*[^/\s.]+[^/]*$#;

my @Samples;
open L,'<',$opt_s or die "[x]Error opening $opt_s: $!\n";
print STDERR "[!]Sample Order: ";
while (<L>) {
	m#([^/]+)/[^/]+$#;
	push @Samples,$1;
	print STDERR (scalar @Samples),':[',$1,"] ";
}
print STDERR "\n";

my ($C_PSNP,$C_iSNP,%SNP)=(0,0);
print STDERR "[!]Parsing SNP ";
open P,'<',$opt_i or die "[x]Error opening $opt_i: $!\n";
while (my $file=<P>) {
	chomp $file;
	open SNP,'<',$file or (warn "\n[!]Error opening $file: $!\n" and next);
	print STDERR ".\b";
	my ($chr,$pos,$ref,$tail,$i);
	while (<SNP>) {
		($chr,$pos,$ref,$tail)=split /\t/;
		++$C_PSNP;
		my @indSNP=split / /,$tail;	# /[ACGTRYMKSWHBVDNX-]/
		for my $s (@Samples) {
			$tail=shift @indSNP;	# ReCycle ...
			if ($tail eq '-') {$tail=undef;}
			 else {++$C_iSNP;}
			$SNP{$chr}{$pos}{$s}=$tail;
		}
	}
	print STDERR '-';
}
warn "[!]PSNP: $C_PSNP, iSNP: $C_iSNP\n";

my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
