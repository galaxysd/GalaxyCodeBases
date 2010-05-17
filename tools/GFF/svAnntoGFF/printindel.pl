#!/usr/bin/perl
use strict;
#use warnings;
use DBI;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.1.0;
$|=1;

our $opts='o:s:d:bv';
our($opt_o, $opt_s, $opt_d, $opt_v, $opt_b);

our $help=<<EOH;
\t-s Specie name for shareing SQLite data file (human)
\t  [Specie name MUST be ths SAME throughout !]
\t-d Result SQLite data file (_result_indel.sqlite)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs

ATTENTION: This script is only for Human !
EOH
# \t-d dbSNP SQLite data file (_tdbSNP.sqlite)

ShowHelp();

$opt_o='__result_indel.txt' if ! $opt_o;
$opt_d='_result_indel.sqlite' if ! $opt_d;
$opt_s='human' if ! $opt_s;

print STDERR "From [$opt_d]\n  to [$opt_o]\tSpecie:[$opt_s]\n";
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
my %attr = (
    RaiseError => 0,
    PrintError => 1,
    AutoCommit => 0
);

my $dbh = DBI->connect('dbi:SQLite:dbname='.$opt_d,'','',\%attr)
 or die $DBI::errstr;

my %CNT;

my $sth = $dbh->prepare( "SELECT chrid,primary_inf,COUNT(1) FROM reindel$opt_s
 GROUP BY chrid,primary_inf ORDER BY primary_inf" );
$sth->execute;

my ($rv,%chrT,%pinfT,$key);
while ($rv=$sth->fetchrow_arrayref) {
	$CNT{$$rv[0]}{$$rv[1]}=$$rv[2];
	$pinfT{SUM}+=$$rv[2];
#	print "[$_] " for @$rv;
#	print "\n";
}

my @chr=(1..22,'X','Y','M');
my @pinf=qw(5-UTR 3-UTR CDS intron UNKNOWN);

for $key (@pinf) {
	$pinfT{$key}+=$CNT{$_}{$key} for @chr;
}
for $key (@chr) {
	$chrT{$key}+=$CNT{$key}{$_} for @pinf;
}

print "\nSummary of InDel\n\nChromosome\tSV_number\t5'-UTR\t3'-UTR\tCDS\tIntron\tnon_Gene\n";
for $key (@chr) {
	print "Chromosome$key\t$chrT{$key}\t";
	for (@pinf) {
		if (defined $CNT{$key}{$_}) {print "$CNT{$key}{$_}\t";}
		 else {print "0\t";}
	}
	print "\n";
}
print "Total     \t";
print "$pinfT{$_}\t" for ('SUM',@pinf);
print "\n";
$dbh->rollback;
$dbh->disconnect;

my $stop_time = [gettimeofday];
#close $infile;
#$|=1;
print STDERR "\n Time Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s)\n";

print STDERR "\033[32;1m Please use [$opt_s] as Specie name in later steps.\033[0;0m\n";
__END__
