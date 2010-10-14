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
\t-i Genotype File List in format: /ChrID, ril file/
\t-o Output File in csvsr format for R/qtl (gen_rot.csv)
\t-v Verbose level (undef=0)
\t-b No pause for batch runs
\t-d Debug Mode, for test only
EOH

ShowHelp();

die "[x]Must specify -i gen_in.lst !\n" unless defined $opt_i;
die "[x]-i $opt_i not exists !\n" unless -f $opt_i;
$opt_o='gen_rot.csv' if ! $opt_o;
no warnings;
$opt_v=int $opt_v;
use warnings;

print STDERR "From [$opt_i] to [$opt_o]\n";
print STDERR "DEBUG Mode on !\n" if $opt_d;
print STDERR "Verbose Mode [$opt_v] !\n" if $opt_v;
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my %TRANS=(
	'0'   => 'A',
	'1'   => 'B',
	'0.5' => 'H',
	'NA'  => '-',
);

my $start_time = [gettimeofday];
#BEGIN
my (@ChrID,%File,@ID);
open RIL,'<',$opt_i  or die "[x]Error opening $opt_i: $!\n";
while (<RIL>) {
	next if /^\s*$/;
	chomp;
	my ($ChrID,$file)=split /\t/;
	push @ChrID,$ChrID;
	$File{$ChrID}=$file;
}
close RIL;
open O,'>',$opt_o  or die "[x]Error opening $opt_o: $!\n";
for my $chr (@ChrID) {
	open I,'<',$File{$chr}  or die "[x]Error opening $File{$chr}: $!\n";
	my $temp='#';
	while ($temp =~ /^#/) {
		$temp=<I>;
	}
	chomp $temp;
	unless (@ID) {
		@ID=split /\s+/,$temp;
		@ID=sort @ID;
		print O 'id,,,',join(',',@ID),"\n";
	} else {
		my @t=split /\s+/,$temp;
		@t=sort @t;
		for (0..$#ID) {
			die "[x]Different Individual ID found in $File{$chr} !\n" if $ID[$_] ne $t[$_];
		}
	}
	while (<I>) {
		next if /^#/;
		chomp;
		my ($Pos,@Dat)=split /\s+/;
		$_=$TRANS{$_} for @Dat;
		print O 'SNP',$chr,'-',$Pos,',',$chr,',0,',join(',',@Dat),"\n";
	}
	close I;
}
close O;
