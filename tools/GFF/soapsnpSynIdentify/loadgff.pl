#!/bin/env perl
use lib '/share/raid010/resequencing/soft/lib';
use lib 'E:/BGI/toGit/perlib/etc';
use strict;
use warnings;
use DBI;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.2.4;

our $opts='i:o:s:bva';
our($opt_i, $opt_o, $opt_s, $opt_v, $opt_b, $opt_a);

our $help=<<EOH;
\t-i Input GFF3 annotation file (human.gff) [ *.gz/*.bz2 is OK ]
\t  [Chromosome name will be striped.]
\t-s Specie of the GFF3 (human)
\t  [Specie name MUST be ths SAME throughout !]
\t-o Output SQLite data file (_dbGFF.sqlite)
\t-a Additive mode to share one SQLite data file among multi-Species
\t   [Output file will be OVERWRITE without -a !]
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_i='human.gff' if ! defined $opt_i;
$opt_o='_dbGFF.sqlite' if ! $opt_o;
$opt_s='human' if ! $opt_s;

print STDERR "From [$opt_i] to [$opt_o], Specie:[$opt_s]\n";
print STDERR "[!]Additive mode SET.\n" if $opt_a;
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];

unlink $opt_o unless $opt_a;

my %attr = (
    RaiseError => 0,
    PrintError => 1,
    AutoCommit => 0
);
my $dbh = DBI->connect('dbi:SQLite:dbname='.$opt_o,'','',\%attr) or die $DBI::errstr;

my $sql=q/
CREATE TABLE IF NOT EXISTS gff{---}
(  chrid TEXT,
   primary_inf TEXT,
   start INTEGER,
   end INTEGER,
   strand TEXT,
   frame TEXT,
   name TEXT,
   groups TEXT );
/;
for (split /;/,$sql) {
	next if /^\s*$/;
	s/{---}/$opt_s/g;
	$dbh->do($_) or die $dbh->errstr;
}
$dbh->commit;

my $sth = $dbh->prepare( "INSERT INTO gff$opt_s ( chrid,primary_inf,start,end,strand,frame,groups,name ) VALUES ( ?,?,?,?,?,?,?,? )" );

my $infile;
if ($opt_i =~ /.bz2$/) {
	open( $infile,"-|","bzip2 -dc $opt_i") or die "Error: $!\n";
} elsif ($opt_i =~ /.gz$/) {
 	open( $infile,"-|","gzip -dc $opt_i") or die "Error: $!\n";
} else {open( $infile,"<",$opt_i) or die "Error: $!\n";}

sub unescape {
  my $v = shift;
  $v =~ tr/+/ /;
  $v =~ s/%([0-9a-fA-F]{2})/chr hex($1)/ge;
  return $v;
}

sub read_gene_gff($$) {
    my $file=$_[0];
    my $dbh=$_[1];
    my @dat;
    while (<$infile>) {
	next if /^\s+$/ or /^;+/;
	chomp;
	s/\r//g;	# There DO be some file with \r. DAMN the MS and APPLE !
	my ($seqname, $source, $primary, $start, $end,
	$score, $strand, $frame, $groups) = split /\t/;	# well, reading file no need to Optimize
	$seqname =~ s/^chr
			(?>
				((?<=^chr)o)?
				((?<=^chro)m)?
				((?<=^chrom)o)?
				((?<=^chromo)s)?
				((?<=^chromos)o)?
				((?<=^chromoso)m)?
				((?<=^chromosom)e)?
			)//xi;
	my @groups = split(/\s*;\s*/, $groups);
	my (%groups,$name);
    for my $group (@groups) {
	my ($tag,$value) = split /=/,$group;
	$tag             = unescape($tag);
	my @values       = map {unescape($_)} split /,/,$value;
	$groups{$tag}=\@values;	# patch for those alter-splices
    }
	my @name_order=qw/Parent ID/;
	@name_order=qw/ID Parent/ if $primary =~ /mRNA/;
	for (@name_order) {
		if ($groups{$_}) {$name=$groups{$_};last;}
	}
	for (@$name) {
		@dat=($seqname,$primary,$start,$end,$strand,$frame,$groups,$_);
		$$dbh->execute( @dat );
		print "$seqname,$primary,$start,$end,$strand,$frame,$_\n" if $opt_v;
	}
    }
    $$dbh->finish;
    return 1;
}

read_gene_gff($infile,\$sth);

close $infile;
$dbh->commit;

$sql=q/
CREATE INDEX IF NOT EXISTS cse{---} ON gff{---}(chrid,start,end);
CREATE INDEX IF NOT EXISTS n{---} ON gff{---}(name);
/;
for (split /;/,$sql) {
	next if /^\s*$/;
	s/{---}/$opt_s/g;
	$dbh->do($_) or warn $dbh->errstr;
}
$dbh->commit;
$dbh->disconnect;

my $read_time = [gettimeofday];
system 'bzip2','-9k',$opt_o;

my $stop_time = [gettimeofday];

$|=1;
print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time ),
	" second(s).\n   Parseing GFF file used\t",tv_interval( $start_time, $read_time ),
	" second(s).\n";

print STDERR "\033[32;1m Please use [$opt_s] as Specie name in later steps.\033[0;0m\n";
