#!/usr/bin/perl -w
use threads;
use strict;
use warnings;
use DBI;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.2.0;

our $opts='i:o:s:bv';
our ($opt_i, $opt_o, $opt_s, $opt_v, $opt_b);

our $help=<<EOH;
\t-i Input SQLite data file (_result.sqlite)
\t-s Specie name (human)
\t  [Specie name MUST be ths SAME throughout !]
\t-o Output txt file (_output_gene.txt)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_o='_output_gene.txt' if ! defined $opt_o;
$opt_i='_result.sqlite' if ! $opt_i;
$opt_s='human' if ! $opt_s;

print STDERR "From [$opt_i] to [$opt_o], Specie:[$opt_s]\n";
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
my %attr = (
    RaiseError => 0,
    PrintError => 1,
    AutoCommit => 0
);
our $dbh = DBI->connect('dbi:SQLite:dbname='.$opt_i,'','',\%attr) or die $DBI::errstr;
my $srdhi = $dbh->prepare("SELECT ALL name,chged FROM res$opt_s
 WHERE primary_inf LIKE '%CDS'") or warn $dbh->errstr;
$srdhi->execute;
###################
	my (%gene_s,%gene_ns,$chged);
	open FH,'>',$opt_o or die "Error: $!\n";
	print FH "#name\tSyn\tNon_Syn\n";
	while (my $ary_ref = $srdhi->fetchrow_arrayref) {
		my ($name,$aa_chg)=@$ary_ref;
		$gene_s{$name}=0 if ! defined $gene_s{$name};
		$gene_ns{$name}=0 if ! defined $gene_ns{$name};
		if ( $aa_chg == 0 ) { $chged=0;++$gene_s{$name} }
		  elsif ( $aa_chg == 1 ) { $chged=1;++$gene_ns{$name} }
print "$name\t$aa_chg\n" if $opt_v;
	}
	for (sort keys %gene_s) {
		print FH "$_\t$gene_s{$_}\t$gene_ns{$_}\n";
	}
	close FH;
