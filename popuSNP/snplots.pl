#!/bin/env perl
use lib '/share/raid010/resequencing/resequencing/tmp/bin/annotation/glfsqlite';
use strict;
use warnings;
use DBI;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.2;

our $opts='i:o:s:f:bv';
our ($opt_i, $opt_o, $opt_s, $opt_v, $opt_b, $opt_l);

our $help=<<EOH;
\t-i PSNP list (./psnp.lst) for chrid.individual.finalSNPs
\t-l merge.list for scaffold positions (./watermelon.merge.list)
\t-s GLF list (./glf.list), will use \$1 of (([^/]+)/[^/]+$) for sample names
\t-o Output Prefix (./indGenomes/ig_)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_o='./indGenomes/ig_' if ! defined $opt_o;
$opt_i='./psnp.lst' if ! $opt_i;
$opt_s='./glf.list' if ! $opt_s;
$opt_l='./watermelon.merge.list' if ! $opt_l;

$opt_i =~ s/\/$//;
print STDERR "From [$opt_i] to [$opt_o], with [$opt_s][$opt_l]/\n";
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

### Begin SQL
my %attr = (
    RaiseError => 0,
    PrintError => 1,
    AutoCommit => 0 );
my $dbh = DBI->connect('dbi:SQLite:dbname=:memory:','','',\%attr) or die $DBI::errstr;

# ChrNew2 scaffold90      1       45293
my $sql=q/
CREATE TABLE IF NOT EXISTS merge
(  chrid TEXT,
   scaffold TEXT,
   start INTEGER,
   end INTEGER  );
/;
# The rowid value can be accessed using one of the special names "ROWID", "OID", or "_ROWID_".
for (split /;/,$sql) {
	next if /^\s*$/;
	$dbh->do($_) or die $dbh->errstr;
}
$dbh->commit;

sub doindex() {
	my $sql=q/
CREATE INDEX IF NOT EXISTS se ON merge(start,end);
CREATE INDEX IF NOT EXISTS c ON merge(chrid);
/;
	for (split /;/,$sql) {
		next if /^\s*$/;
		$dbh->do($_) or warn $dbh->errstr;
	}
	$dbh->commit;
}

my $sthi = $dbh->prepare( 'INSERT INTO merge ( chrid,scaffold,start,end ) VALUES ( ?,?,?,? )' );
my $stho = $dbh->prepare( 'SELECT DISTINCT scaffold,start,end FROM merge WHERE chrid=? AND ? BETWEEN start AND end' );
#$|=1;

my @Scaffolds;
open L,'<',$opt_l or die "[x]Error opening $opt_l: $!\n";
print STDERR "[!]Scaffolds: ";
my $i=0;
while (<L>) {
	chomp;
	my ($chrid,$scaffold,$start,$end)=split /\t/;
	$sthi->execute($chrid,$scaffold,$start,$end);
	++$i;
}
close L;
&doindex;
print STDERR "$i loaded.\n";

print STDERR "[!]Parsing SNP ";
open P,'<',$opt_i or die "[x]Error opening $opt_i: $!\n";



my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
