#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BIOPIC <galaxy001@gmail.com>
Version: 1.0.0 @ 20120720
=cut
use strict;
use warnings;
use Vcf;
#use GTF;
use DBI;
use Galaxy::IO;
use Galaxy::ChromString;
use Data::Dump qw(ddx);

die "Usage: $0 <db_file> <bwa ann file> <infile> <out>\n" if @ARGV<4;
my $dbfile = shift;
my $annfile = shift;
my $infs = shift;
my $outfs = shift;

warn "From [$infs] to [$outfs] with [$dbfile,$annfile]\n";

open D,'<',$dbfile or die;
open A,'<',$annfile or die;
my $IN = openfile($infs);
open O,'>',$outfs or die;
print O "# From [$infs] to [$outfs] with [$dbfile,$annfile]\n";

my $dbh = DBI->connect(
    "dbi:SQLite:dbname=mem.sqlite","","", 
    {RaiseError => 0, PrintError => 1, AutoCommit => 0}
) or die $DBI::errstr;
my $sql=q/
CREATE TABLE IF NOT EXISTS gff
(  chrid TEXT,
   primary_inf TEXT,
   start INTEGER,
   end INTEGER,
   strand TEXT,
   frame TEXT,
   name TEXT );
/;
for (split /;/,$sql) {
        next if /^\s*$/;
        $dbh->do($_) or die $dbh->errstr;
}
$dbh->commit;

my $sth = $dbh->prepare( "INSERT INTO gff ( chrid,primary_inf,start,end,strand,frame,name ) VALUES ( ?,?,?,?,?,?,? )" );


my (%ChrLen,$ChrDat,@IDs,@Lens);
<A>;
while (<A>) {
	my (undef,$id) = split / /;	# gi|477502308|ref|NW_004457745.1|
	if ($id =~ /ref\|([^|]+)\|/) {
		$id = $1;
	}
	$_ = <A>;
	my (undef,$len) = split / /;
	$ChrLen{$id} = $len;
	push @IDs,$id;
	push @Lens,$len;
	print "$id, $len\n";
}
$ChrDat = ChrStrInit(\@IDs,\@Lens);
close A;

my %bitflag = (
		'CDS' => 1,
		'mRNA' => 2
	);

my $secname = ']';
my (%GeneDat, %Gene2Chr, $strand, $seqname, $flag);
while (<D>) {
	next if /^(#|((\s)*$))/;
	my ($primary, $start, $end, $frame);
	if (/^\[([^]]*)\] ([+-]) ([^ ]+) (.+)$/) {
#ddx $Gene2Chr{$secname};
		$secname = $1;
		$strand = $2;
		$seqname = $3;	# NW_004465862.1
		$flag = $4;
		#print "[$secname, $strand, $seqname, $flag]\n";
		die "[$_]" if length $secname == 0;
	} else {
		chomp;
		($primary, $start, $end, $frame) = split /\t/,$_;
		if ( $primary eq 'CDS' or $primary eq 'mRNA' ) {
			push @{$GeneDat{"$secname\n$strand"}{$primary}},[$start, $end, $frame];
			++$Gene2Chr{$secname}{$seqname};
			$sth->execute( $seqname, $primary, $start, $end, $strand, $frame, $secname );
=pod
			for my $i ( $start .. $end ) {
				if ( ChrStrANDBase($ChrDat, $seqname, $i, $bitflag{$primary}) > 0 ) {
					ChrStrORBase($ChrDat, $seqname, $i, ($bitflag{$primary} << 4) );
				}
				ChrStrORBase($ChrDat, $seqname, $i, $bitflag{$primary});
			}
=cut
		}
#ddx $GeneDat{"$secname\n$strand"};
	}
	#ddx $$ChrDat{$seqname};
}
close D;
$dbh->commit;

$sql=q/
CREATE INDEX IF NOT EXISTS cseGFF ON gff(chrid,start,end);
CREATE INDEX IF NOT EXISTS nGFF ON gff(name);
/;
for (split /;/,$sql) {
        next if /^\s*$/;
        $dbh->do($_) or warn $dbh->errstr;
}
$dbh->commit;
my $sthi = $dbh->prepare( "SELECT * FROM gff WHERE chrid=? AND ? BETWEEN start AND end" );

my @bamnames;
while (<$IN>) {
	if (/^#/) {
		if (/^#CHROM\tPOS\t/) {
			chomp;
			@bamnames = split /\t/;
			@bamnames = splice @bamnames, 9;
print "[@bamnames]\n";
		}
		next;
	}
	chomp;
	my ($chr, $pos, undef, $refbase, $altbases, $QUAL, undef, $INFO, $FORMAT, @data) = split /\t/;
}
close $IN;

close O;

__END__
perl ano_zones.pl Dasnov3.db Dasnov3.ann bam2.bcgv.vcf.gz t
