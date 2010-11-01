#!/bin/env perl
use lib '/share/raid010/resequencing/soft/lib';
use lib 'E:/BGI/toGit/perlib/etc';
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
#use File::Basename;
use Galaxy::ShowHelp;
#use FindBin qw($RealBin);

$main::VERSION=0.0.1;
#my $SCRIPTS="$RealBin/../scripts";

our $opts='i:o:v:bd';
our($opt_i, $opt_o, $opt_v, $opt_b, $opt_d);

#our $desc='';
our $help=<<EOH;
\t-i Phenotype Data File in format : /id, phe1, phe2, .../ with Title
\t-o Output File in csvsr format for R/qtl (phe_rot.csv)
\t-v Verbose level (undef=0)
\t-b No pause for batch runs
\t-d Debug Mode, for test only
EOH

ShowHelp();

die "[x]Must specify -i phe_in.csv !\n" unless defined $opt_i;
die "[x]-i $opt_i not exists !\n" unless -f $opt_i;
$opt_o='phe_rot.csv' if ! $opt_o;
no warnings;
$opt_v=int $opt_v;
use warnings;

print STDERR "From [$opt_i] to [$opt_o]\n";
print STDERR "DEBUG Mode on !\n" if $opt_d;
print STDERR "Verbose Mode [$opt_v] !\n" if $opt_v;
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
my (@ID,%Phe);
open I,'<',$opt_i  or die "[x]Error opening $opt_i: $!\n";
my $line=<I>;
$line =~ s/ *\r//g;
chomp $line;
my @phe=split /\s*\,\s*/,$line;
shift @phe;
while (<I>) {
	next if /^;/;
	s/ *\r//g;
	chomp;
	s/\,\s*$//;
	my @dat=split /\s*\,\s*/;
	my $id = shift @dat;
	push @ID,$id;
	$Phe{$id}=\@dat;
}
@ID = sort @ID;
close I;
open O,'>',$opt_o  or die "[x]Error opening $opt_o: $!\n";
for my $i (0..$#phe) {
	print O $phe[$i];
	for my $id (@ID) {
		print O ',',$Phe{$id}->[$i];
	}
	print O "\n";
}
print O 'id,',join(',',@ID),"\n";
close O;
