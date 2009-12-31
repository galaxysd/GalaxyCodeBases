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
\t-i Genome SQLite data file (_dbGFF.sqlite)
\t-s Specie name (human)
\t  [Specie name MUST be ths SAME throughout !]
\t-o Output log file (_dbGFF.check) [not Working now]
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_o='_dbGFF.check' if ! $opt_o;
$opt_i='_dbGFF.sqlite' if ! $opt_i;
$opt_s='human' if ! $opt_s;

print STDERR "From [$opt_i] to [$opt_o], Specie:[$opt_s]\n";
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];

my %attr = (
    RaiseError => 0,
    PrintError => 1,
    AutoCommit => 0
);
my $dbh = DBI->connect('dbi:SQLite:dbname='.$opt_i,'','',\%attr) or die $DBI::errstr;

my $sth = $dbh->prepare( "SELECT name,strand FROM gff$opt_s
 WHERE primary_inf LIKE '%mRNA' ORDER BY chrid,start,strand" );
my $sthp = $dbh->prepare( "SELECT chrid,start,end,strand,frame FROM gff$opt_s
 WHERE name=? AND primary_inf LIKE '%CDS' ORDER BY start ASC" );
my $sthm = $dbh->prepare( "SELECT chrid,start,end,strand,frame FROM gff$opt_s
 WHERE name=? AND primary_inf LIKE '%CDS' ORDER BY start DESC" );
$sth->execute;
my ($rv,$res,%COUNT,$chrid,$frameC);	#$last_b,$last_e,$strand,$frame,$last_f,
$COUNT{'1gene'}=$COUNT{'2right'}=$COUNT{'3opp'}=$COUNT{'4err'}=0;
while ($rv=$sth->fetchrow_arrayref) {
	++$COUNT{'1gene'};
print @$rv,"\n" if $opt_v;
	if ($$rv[1] eq '+') {$sthp->execute($$rv[0]);$res=$sthp->fetchall_arrayref;}
	  elsif ($$rv[1] eq '-') {$sthm->execute($$rv[0]);$res=$sthm->fetchall_arrayref;}
	  else {warn "Unexpected strand: $$rv[1] !\n"; next;}
	my $last_res=shift @$res;	# chrid,start,end,strand,frame
	next unless defined $last_res;
	my $last_f=0;
print "[!]@$rv:\t@$last_res\t",1+$$last_res[2]-$$last_res[1],"\n" if $$last_res[4] ne '0';
	++$COUNT{'4err'} if $$last_res[4] != 0;
print "@$last_res\t",1+$$last_res[2]-$$last_res[1],"\n" if $opt_v;
	for (@$res) {
		$last_f=$frameC=(3-((1+abs($$last_res[1]-$$last_res[2])-$last_f)%3))%3;
print "@$_\t",1+$$_[2]-$$_[1],"\t$frameC\n" if $opt_v;
		my $frame=$$_[4];
		$last_res=$_;
		next if $frameC == 0 and $frame==0;
		if ($frameC == $frame) {++$COUNT{'2right'};}
		  elsif (abs($frameC-$frame)==1) {++$COUNT{'3opp'};}
		  else {++$COUNT{'4err'};}
	}
}
for (sort keys %COUNT) {
#	$COUNT{$_}=0 unless $COUNT{$_};
	print "$_:\t$COUNT{$_}\n"
}

if ($COUNT{'4err'} > 0) {print "Wrong GFF frames !\nPlease run fixgff.pl\n"}
  elsif ($COUNT{'2right'} > 0 and  $COUNT{'3opp'} == 0) {print "Stranded GFF frames ! No fix needed.\n"}
  elsif ($COUNT{'2right'} == 0 and  $COUNT{'3opp'} > 0) {print "BGI GFF frames !\nPlease run fixgff.pl or you have to use -f in update_aa.pl\n"}
  else {print "Mixed GFF frames.\nPlease run fixgff.pl\n"}
$dbh->rollback;
$dbh->disconnect;

my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";

print STDERR "\033[32;1m Please use [$opt_s] as Specie name in later steps.\033[0;0m\n";
