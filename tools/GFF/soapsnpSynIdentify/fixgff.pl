#!/bin/env perl
use lib '/share/raid010/resequencing/soft/lib';
use lib 'E:/BGI/toGit/perlib/etc';
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
\t-o Output SQLite data file (_dbGFF.fixed.sqlite)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_o='_dbGFF.fixed.sqlite' if ! $opt_o;
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

system 'cp','-f',$opt_i,$opt_o and die "[x]Cannot initialize output file $opt_o !\n";;
my $dbh = DBI->connect('dbi:SQLite:dbname='.$opt_o,'','',\%attr) or die $DBI::errstr;

my $sth = $dbh->prepare( "SELECT name,strand FROM gff$opt_s
 WHERE primary_inf LIKE '%mRNA' ORDER BY chrid,start,strand" );
my $sthp = $dbh->prepare( "SELECT chrid,start,end,strand,frame FROM gff$opt_s
 WHERE name=? AND primary_inf LIKE '%CDS' ORDER BY start ASC" );
my $sthm = $dbh->prepare( "SELECT chrid,start,end,strand,frame FROM gff$opt_s
 WHERE name=? AND primary_inf LIKE '%CDS' ORDER BY start DESC" );
$sth->execute;

my $sthf = $dbh->prepare( "UPDATE gff$opt_s SET frame=?
 WHERE name=? AND chrid=? AND start=? AND end=? AND primary_inf LIKE '%CDS'" );

my ($rv,$res,%COUNT,$chrid,$frameC);	#$last_b,$last_e,$strand,$frame,$last_f,
$COUNT{'1gene'}=$COUNT{'2right'}=$COUNT{'3opp'}=$COUNT{'4err'}=0;
while ($rv=$sth->fetchrow_arrayref) {
	++$COUNT{'1gene'};
	if ($$rv[1] eq '+') {$sthp->execute($$rv[0]);$res=$sthp->fetchall_arrayref;}
	  elsif ($$rv[1] eq '-') {$sthm->execute($$rv[0]);$res=$sthm->fetchall_arrayref;}
	  else {warn "Unexpected strand: $$rv[1] !\n"; next;}
	my $last_res=shift @$res;	# chrid,start,end,strand,frame
	next unless defined $last_res;
	my $last_f=0;
	$sthf->execute(0,$$rv[0],$$last_res[0],$$last_res[1],$$last_res[2]) if $$last_res[4] ne '0';	# some gff come with '.'
	# first is always 0. There DO be such wrong file.
	for (@$res) {
		$last_f=$frameC=(3-((1+abs($$last_res[1]-$$last_res[2])-$last_f)%3))%3;
		my $frame=$$_[4];
		$last_res=$_;
		$sthf->execute($frameC,$$rv[0],$$_[0],$$_[1],$$_[2]);
		next if $frameC == 0;
		unless ($frame =~ /\d+/) { ++$COUNT{'5undef'};next; }   # fix for Frame eq '.'
		if ($frameC == $frame) {++$COUNT{'2right'};}
		  elsif (abs($frameC-$frame)==1) {++$COUNT{'3opp'};}
		  else {++$COUNT{'4err'};}
	}
}
for (sort keys %COUNT) {
	print "$_:\t$COUNT{$_}\n"
}

if ($COUNT{'4err'} > 0) {print "Error GFF frames !\n"}
  elsif ($COUNT{'2right'} > 0 and  $COUNT{'3opp'} == 0) {print "Stranded GFF frames !\n"}
  elsif ($COUNT{'2right'} == 0 and  $COUNT{'3opp'} > 0) {print "BGI GFF frames !\n"}
  else {print "???\n"}
print "\nFile fixed !\n";
$dbh->commit;
$dbh->disconnect;

system 'bzip2','-9k',$opt_o;

#END
my $stop_time = [gettimeofday];

$|=1;
print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time ),
	" second(s).\n";

print STDERR "\033[32;1m Please use [$opt_s] as Specie name in later steps.\033[0;0m\n";
__END__

./update_aa.pl -b -o _result1.sqlite > _result1.txt 2>_result1.err &
