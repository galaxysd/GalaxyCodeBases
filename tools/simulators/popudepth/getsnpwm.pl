#!/usr/bin/perl
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.1.1;

our $opts='i:o:bv';
our($opt_i, $opt_o, $opt_v,$opt_b);

our $help=<<EOH;
\t-i Add_ref files list (./indsnp.lst)
\t-o Output PWM file (snp.pwm)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH
# \t-d dbSNP SQLite data file (_tdbSNP.sqlite)

ShowHelp();

$opt_i='./indsnp.lst' if ! defined $opt_i;
$opt_o='snp.pwm' if ! $opt_o;

print STDERR "From [$opt_i] to [$opt_o]\n";
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

#http://doc.bioperl.org/bioperl-live/Bio/Tools/IUPAC.html#BEGIN1
my %IUB = ( A => [qw(A)],
	     C => [qw(C)],
	     G => [qw(G)],
	     T => [qw(T)],
	     U => [qw(U)],
	     M => [qw(A C)],
	     R => [qw(A G)],
	     W => [qw(A T)],
	     S => [qw(C G)],
	     Y => [qw(C T)],
	     K => [qw(G T)],
	     V => [qw(A C G)],
	     H => [qw(A C T)],
	     D => [qw(A G T)],
	     B => [qw(C G T)],
	     X => [qw(G A T C)],
	     N => [qw(G A T C)]
	     );

my $start_time = [gettimeofday];

my @files;
open L,'<',$opt_i or die "Error opening $opt_i: $!\n";
while (<L>) {
	chomp;
	if (-s $_) {push @files,$_;}
	 else {warn "[!]$_ not available.\n";}
}
close L;
warn '[!]Total: ',@files+1," files.\n";

my $snpcount;
open O,'>',$opt_o or die "Error opening $opt_o: $!\n";
for my $file (@files) {
	open I,'<',$opt_i or die "Error opening $opt_i: $!\n";
	while (<I>) {
		chomp;
		my ($chr,$pos,$bases) = split /\t/;
		++$snpcount;
		print O "$chr\t$pos";
		my @base = split / /,$bases;
		my ($i,%counter);
		for (@base) {
			my $arr=$IUB{$_} or next;
			for (@$arr) {
				++$counter{$_};
				++$i;
			}
		}
		for (sort {$counter{$b} <=> $counter{$a}} keys %counter) {
			print O "\t$_:",$counter{$_}/$i;
		}
		print O "\n";
	}
	close I;
}
close O;

warn "[!]Total: $snpcount SNPs.\n[!]Done !\n";
