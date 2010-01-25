#!/bin/env perl
use lib '/share/raid010/resequencing/resequencing/tmp/bin/annotation/glfsqlite';
use lib 'E:\BGI\toGit\perlib\etc';
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;
use GalaxyXS::ChromByte 1.02;
use DBI;

$main::VERSION=0.0.3;

our $opts='i:o:l:t:c:sbv';
our ($opt_i, $opt_o, $opt_l, $opt_v, $opt_b, $opt_s, $opt_t, $opt_c);

our $help=<<EOH;
\t-i Input Genome sequence file (human.fa)
\t-l kmer length (33)
\t-o Output Stat (stat.txt)
\t-s Sort Output by Count DESCending [will be slower]
\t-c Cache size [GB in Memory] (4)
\t-t tmpfile for swap, will be removed (./._tmp_)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_i='human.fa' if ! defined $opt_i;
$opt_o='stat.txt' if ! $opt_o;
$opt_t='./._tmp_' if ! $opt_t;
$opt_c=4 if ! $opt_c;
$opt_c=0.002 if $opt_c<0.002;
$opt_l=33 if ! $opt_l;
$opt_l=int($opt_l);
$opt_l=1 if $opt_l<1;

print STDERR "From [$opt_i] with [$opt_l][$opt_t][$opt_c]G to [$opt_o]\n";
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

open( T,'>',$opt_t) or die "[x]Error creating swap file: $!\n";
close T;
unlink $opt_t;
my $start_time = [gettimeofday];

my %attr = (
    RaiseError => 0,
    PrintError => 1,
    AutoCommit => 0,
);
my $dbh = DBI->connect('dbi:SQLite:dbname='.$opt_t,undef,undef,\%attr) or die $DBI::errstr;
my $sql=q/
PRAGMA synchronous = OFF;
PRAGMA journal_mode = OFF;
PRAGMA temp_store = MEMORY;
PRAGMA auto_vacuum = NONE;
PRAGMA cache_size = /.int(1000000*$opt_c).	# 2000 for 2MB. Bigger, better ?
q/;
CREATE TABLE kmers ( seq TEXT );
/;
for (split /;/,$sql) {
	next if /^\s*$/;
	$dbh->do($_) or die $dbh->errstr,": $_";
}
my $sth = $dbh->prepare( 'INSERT INTO kmers ( seq ) VALUES ( ? );' );

open( IN,'<',$opt_i) or die "[x]Error: $!\n";
my ($total,$seq,$len)=(0);
while (<IN>) {
	s/^>//;
	my $title = $_;
	my $seqname = $1 if($title =~ /^(\S+)/);
	print STDERR "loading >$seqname\t";
	$/=">";
	$seq=<IN>;
	chomp $seq;
	$seq=~s/\s//g;
	$/="\n";
	# Begin Kmers
	$len=length $seq;
	$seq=uc $seq;
#$len=10000000+$opt_l-1;
	print STDERR " :$len -> ";
	#next if $len < $opt_l;	# We already while it
	my $i;
	my $t=$len-$opt_l;	# Will perl do this automatically?
	for ($i=0;$i<=$t;$i++) {
		my $str=substr $seq,$i,$opt_l;
#die if length $str != $opt_l;	# just test
		$sth->execute($str);# or warn $sth->errstr;
	}
	print STDERR $i,"\n";
	$total += $i;
	#$dbh->commit;
	# End Kmers
}
close IN;
my $read_time = [gettimeofday];
$seq=undef;
warn "[!]$total Kmers Inserted !\n";

$dbh->do('CREATE INDEX c ON kmers(seq);') or warn $dbh->errstr,": $_";	# without INDEX, much more memory required for SELECT
my $i_time = [gettimeofday];

if ($opt_s) { $sql='SELECT seq,COUNT(seq) AS count FROM kmers GROUP BY seq ORDER BY count DESC;'; }
 else { $sql='SELECT seq,COUNT(seq) AS count FROM kmers GROUP BY seq;'; }
my $gth = $dbh->prepare( $sql );
$gth->execute;
my $work_time = [gettimeofday];
warn "[!]Count done !\n";

open( O,'>',$opt_o) or die "[x]Error creating output file: $!\n";
while (my $ary_ref = $gth->fetchrow_arrayref) {
	#my ($seq,$count)=@$ary_ref;
	print O $$ary_ref[0],"\t",$$ary_ref[1],"\n";
}
close O;
my $end_time = [gettimeofday];
warn "[!]Output done !\n";

#$dbh->commit;
$dbh->disconnect;
unlink $opt_t;
my $stop_time = [gettimeofday];

$|=1;
print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time ),
	" second(s).\n   Insert used:\t",tv_interval( $start_time, $read_time ),
	" second(s).\n   Index used:\t",tv_interval( $read_time, $i_time ),
	" second(s).\n   Count used:\t",tv_interval( $i_time, $work_time ),
	" second(s).\n   Output used:\t",tv_interval( $work_time, $end_time ),
	" second(s).\n   Clean used:\t",tv_interval( $end_time, $stop_time )," second(s).\n";
__END__

10000000*2 = 20 M

with sort:
Time Elapsed:   653.591479 second(s).
   INSERT used: 190.030188 second(s).
   Count used:  325.451145 second(s).
   Output used: 137.574869 second(s).
   Clean used:  0.535277 second(s).

without sort:
Time Elapsed:   544.888772 second(s).
   INSERT used: 189.922403 second(s).
   Count used:  199.940926 second(s).
   Output used: 154.487606 second(s).
   Clean used:  0.537837 second(s).

without sort, PRAGMA journal_mode = MEMORY;:
Time Elapsed:   547.098508 second(s).
   Insert used: 190.366605 second(s).
   Count used:  200.026188 second(s).
   Output used: 156.246296 second(s).
   Clean used:  0.459419 second(s).

without sort, cache=1 GB:
Time Elapsed:   547.102521 second(s).
   Insert used: 189.94249 second(s).
   Count used:  200.629807 second(s).
   Output used: 155.977248 second(s).
   Clean used:  0.552976 second(s).

2 M:
without sort, cache=10 MB:
Time Elapsed:   56.378994 second(s).
   Insert used: 23.819828 second(s).
   Count used:  16.546491 second(s).
   Output used: 15.898938 second(s).
   Clean used:  0.113737 second(s).

without sort, cache=100 MB:
Time Elapsed:   162.321573 second(s).
   Insert used: 23.709979 second(s).
   Count used:  102.493035 second(s).
   Output used: 36.062866 second(s).
   Clean used:  0.055693 second(s).
