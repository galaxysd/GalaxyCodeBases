#!/bin/env perl
#use threads;
use strict;
use warnings;
use DBI;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;
$main::VERSION=0.1.1;

our $opts='i:o:m:bvp';
our($opt_i, $opt_o, $opt_v, $opt_b, $opt_m, $opt_p);

#our $desc='SoapSort library PCR PE Duplication Remover & Merger (Atom Edition)';
our $help=<<EOH;
\t-i Input file, in format: /^ChrID\\tPos\\t.*\\n\$/
\t-m merge list, in format: /^ChrID\\tscaffoldID\\tChrStart\\tChrStop\\n\$/
\t-o Output file
\t-p only part of the ChrID needs to change
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

warn "[x]Must specify -i !\n" unless $opt_i;
warn "[x]Must specify -m !\n" unless $opt_m;
warn "[x]Must specify -o !\n" unless $opt_o;
exit 1 unless $opt_i and $opt_o and $opt_m;

print STDERR "From [$opt_i] with [$opt_m] to [$opt_o]\n";
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
unless (-s $opt_i) {die "[x]Soaplist [$opt_i] is nothing !\n";}

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
$|=1;

open SAMPLE,'-|',"gzip -dc $opt_m" or die "Error opening $opt_m: $!\n";
while (<SAMPLE>) {
	chomp;
	my ($chrid,$scaffold,$start,$end)=split /\t/;
	$sthi->execute($chrid,$scaffold,$start,$end);
}
close SAMPLE;
&doindex;

open IN,'<',$opt_i or die "Error opening $opt_i: $!\n";
open OUT,'>',$opt_o or die "Error opening $opt_o: $!\n";
my ($chrid,$pos,$scaffold,$start,$end,$qres,$newpos,@last);
while (<IN>) {
# ChrNew1	121900	T T T T W W T T T T T T T T T T W W T T T T T T T - T T T T T T T T T T
	($chrid,$pos,@last)=split /\t/;
	$stho->execute($chrid,$pos);
	$qres = $stho->fetchall_arrayref;
	if ($#$qres == -1) {
		if ($opt_p) {
			print OUT join("\t",$chrid,$pos,@last);
			next;
		} else {
			warn "No info. for $chrid:$pos, skipped.\n";
			next;
		}
	} elsif ($#$qres != 0) { warn "More than 1 hit, using first.\n" }
	($scaffold,$start,$end)=@{$$qres[0]};
	$newpos=1+$pos-$start;
	print OUT join("\t",$scaffold,$newpos,@last);
}
close OUT;
close IN;
#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__

