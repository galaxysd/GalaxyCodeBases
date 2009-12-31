#!/bin/env perl
use strict;
use warnings;
use DBI;
use integer;

unless (@ARGV){
	print "perl $0 <GFF in> <GFF out>\n";
	exit;
}

my $in = shift;
my $out = shift;

### Begin SQL
my %attr = (
    RaiseError => 0,
    PrintError => 1,
    AutoCommit => 0 );
my $dbh = DBI->connect('dbi:SQLite:dbname=:memory:','','',\%attr) or die $DBI::errstr;

my $sql=q/
CREATE TABLE IF NOT EXISTS gff
(  seqname TEXT,
   source TEXT,
   primary_inf TEXT,
   start INTEGER,
   end INTEGER,
   score TEXT,
   strand TEXT,
   frame TEXT,
   groups TEXT );
/;
for (split /;/,$sql) {
	next if /^\s*$/;
	$dbh->do($_) or die $dbh->errstr;
}
$dbh->commit;

my $sthi = $dbh->prepare( 'INSERT INTO gff ( seqname,source,primary_inf,start,end,score,strand,frame,groups ) VALUES ( ?,?,?,?,?,?,?,?,? )' );
my $sthc = $dbh->prepare( 'SELECT DISTINCT seqname,start,end FROM gff WHERE primary_inf = ? ORDER BY start' );
my $stho = $dbh->prepare( 'SELECT DISTINCT primary_inf,start,end FROM gff WHERE seqname=$3 AND (NOT primary_inf=$4) AND (start BETWEEN $1 AND $2) AND (end BETWEEN $1 AND $2)' );

sub Subtraction($$) {
	my ($source,$tosub)=@_;
	my %tmp;
	++$tmp{$_} for @$source;
	delete $tmp{$_} for @$tosub;
	my @tmp=keys %tmp;
	return \@tmp;
}

sub Zones($$$) {
	my ($start,$end,$subzones)=@_;
	my $len=$end-$start;
}

my %primarys;
open IN, "$in" || die "$!\n";
while (<IN>) {
	chomp;
	my ($seqname, $source, $primary, $start, $end,
	$score, $strand, $frame, $groups) = split /\t/;
	++$primarys{$primary};
	$sthi->execute($seqname, $source, $primary, $start, $end,$score, $strand, $frame, $groups);
}
$sql=q/
CREATE INDEX IF NOT EXISTS spepos ON gff(seqname,primary_inf,start,end);
CREATE INDEX IF NOT EXISTS pis ON gff(primary_inf,start);
/;
for (split /;/,$sql) {
	next if /^\s*$/;
	$dbh->do($_) or warn $dbh->errstr;
	}
$dbh->commit;

my %BITs= (
	UTR5	=> 1,
	UTR3	=> 2,
	UTRx	=> 3,	# 1+2
	CDS	=> 4,
	Exon	=> 7,	# 1+2+4
	Intron	=> 8,
	mRNA	=> 15,	# 1+2+4+8
	Promoter	=> 16,
	Gene	=> 9223372036854775807,	# 7FFF FFFF FFFF FFFF = signed long
);	# 63-5=58 bits for other attr

my @KEYs=keys %primarys;
my @CDS=grep /cds/i,@KEYs;
my @Exon=grep /exon/i,@KEYs;
my @Intron=grep /intron/i,@KEYs;
my @mRNA=grep /mRNA/i,@KEYs;
my @Gene=grep /gene/i,@KEYs;
my @UTR=grep /UTR/i,@KEYs;
my @UTR5=grep /5|five/i,@UTR;
my @UTR3=grep /3|three/i,@UTR;
my @UTRx=@{&Subtraction(\@UTR,[@UTR5,@UTR3])};
=pod
[five_prime_utr mRNA gene three_prime_utr exon CDS UTR]
=cut
#print "[@keys]\n[@CDS][@exon][@intron][@mRNA][@gene][@UTR][@UTR5][@UTR3][@UTRx]\n";

seek IN,0,0;
my $a=<IN>;
print $a;