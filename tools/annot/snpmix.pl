#!/bin/env perl
use lib '/share/raid010/resequencing/resequencing/tmp/bin/annotation/glfsqlite';
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.1;

our $opts='i:o:s:bv:f:';
our ($opt_i, $opt_o, $opt_s, $opt_v, $opt_b, $opt_f);

our $help=<<EOH;
\t-i PSNP list (./psnp.lst) for chrid.individual.finalSNPs
\t-f fabyChr path (./faByChr/) for [chr].fa
\t-s GLF list (./glf.list), will use \$1 of /([^/]+)/[^/]+$/ for sample names
\t-o Output Prefix (./indGenomes/ig_)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_o='./indGenomes/ig_' if ! defined $opt_o;
$opt_i='./psnp.lst' if ! $opt_i;
$opt_s='./glf.list' if ! $opt_s;
$opt_f='./faByChr/' if ! $opt_f;

$opt_i =~ s/\/$//;
$opt_f =~ s/\/$//;
print STDERR "From [$opt_i] to [$opt_o], with [$opt_s][$opt_f]/\n";
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];

system('mkdir','-p',$opt_o);
system('rmdir',$opt_o) if $opt_o =~ #/[^/\s.]+[^/]*$#;

my @Samples;
open L,'<',$opt_s or die "Error opening $opt_s: $!\n";
print STDERR "[!]Sample Order: ";
while (<L>) {
	/([^/]+)/[^/]+$/;
	push @Samples,$1;
	print STDERR $1,"\t";
}
print STDERR "\n";

open P,'<',$opt_i or die "Error opening $opt_i: $!\n";
while my $file (<P>) {
	chomp $file;
	my %SNP;
	open SNP,'<',$file or (warn "Error opening $file: $!\n" and next);
	while (SNP) {
		my ($chr,$pos,$ref,$tail)=split /\t/;
		my @indSNP;
		for (split / /,$tail) {
			next unless /[ACGTRYMKSWHBVDNX-]/;
			push @indSNP,$_;
		}
		$SNP{$pos}=[$ref,\@indSNP];
	}
}



my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
