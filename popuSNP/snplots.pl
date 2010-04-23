#!/bin/env perl
use lib '/share/raid010/resequencing/resequencing/tmp/bin/annotation/glfsqlite';
use strict;
use warnings;
use DBI;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.2;

our $opts='i:o:s:l:m:bv';
our ($opt_i, $opt_o, $opt_s, $opt_v, $opt_b, $opt_l, $opt_m);

our $help=<<EOH;
\t-i PSNP list (./npsnp.lst) for chrid.individual.finalSNPs
\t-l scaffold.list (./scaffold.N90.id)
\t-m merge.list for scaffold positions (./watermelon.merge.list)
\t-s GLF list (./glf.list), will use \$1 of (([^/]+)/[^/]+$) for sample names
\t-o Output Prefix (./plots/p_)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_o='./plots/p_' unless defined $opt_o;
$opt_i='./npsnp.lst' unless $opt_i;
$opt_s='./glf.list' unless $opt_s;
$opt_l='./scaffold.N90.id' unless $opt_l;
$opt_m='./watermelon.merge.list' unless $opt_m;

$opt_i =~ s/\/$//;
print STDERR "From [$opt_i] to [$opt_o], with [$opt_s][$opt_l][$opt_m]\n";
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
close L;
print STDERR "\n";

my ($LEN,%Scaffolds,@Scaffoldo)=(0);
open L,'<',$opt_l or die "[x]Error opening $opt_l: $!\n";
print STDERR "[!]Scaffolds: ";
while (<L>) {
	chomp;
	--$Scaffolds{$_};
	#print STDERR (scalar (keys %Scaffolds)),':[',$_,"] ";
}
close L;
print STDERR scalar (keys %Scaffolds),"\n";

open L,'<',$opt_m or die "[x]Error opening $opt_m: $!\n";
print STDERR "[!]Scaffolds Len: ";
my ($i,$len)=(0);
while (<L>) {
	chomp;
	my ($chrid,$scaffold,$start,$end)=split /\t/;
	next unless exists $Scaffolds{$scaffold};
	$len=$end-$start+1;
	$Scaffolds{$scaffold}=[$len,0];
	$LEN += $len;
warn "$chrid,$scaffold,$start,$end\t$len\t$LEN\n" if $opt_v;
	++$i;
}
close L;
print STDERR "$LEN\n";
@Scaffoldo = sort {${$Scaffolds{$b}}[0] <=> ${$Scaffolds{$a}}[0]} keys %Scaffolds;
my $i=0;
for (@Scaffoldo) {
	${$Scaffolds{$_}}[1]=$i;
	$i += ${$Scaffolds{$_}}[0];
warn "$_\t${$Scaffolds{$_}}[0]\t${$Scaffolds{$_}}[1]\n"  if $opt_v;
}

my (%Plot,$POS);
print STDERR "[!]Parsing SNP ";
open P,'<',$opt_i or die "[x]Error opening $opt_i: $!\n";
while (my $file=<P>) {
	chomp $file;
	open SNP,'<',$file or (warn "\n[!]Error opening $file: $!\n" and next);
	print STDERR ".\b";
	while (<SNP>) {
		my ($scaffold,$pos,$ref,$tail)=split /\t/;
		next unless exists $Scaffolds{$scaffold};
		$POS = ${$Scaffolds{$_}}[1] + $pos;
		my @Values;
		for (split / /,$tail) {
			next unless /[ACGTRYMKSWHBVDNX-]/;
			if (/[ACGTRYMKSWHBVDNX]/) {
				push @Values,;
			} elsif ($_ eq '-') {
				;
			} else {next;}
		}
	}
	close SNP;
}
close P;

my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";

