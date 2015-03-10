#!/usr/bin/env perl
use strict;
use warnings;
use DBI;
use Data::Dump qw(ddx);

# http://www.perlmonks.org/?node_id=393426
use FindBin;
use lib $FindBin::Bin;
use extractfuse;

die "Usage: $0 <ins_size> <merged.sqlite> <out_prefix>\n" if @ARGV < 3;
my ($isize,$in,$out)=@ARGV;

our $DEBUG=2.9;
our $minLen = 5;

our (%Genome,%ChrLen);	# Let's share them.

my $dbh = DBI->connect("dbi:SQLite:dbname=./$in","","",{
	RaiseError => 1,
	PrintError => 1,
	AutoCommit => 0
}) or die $DBI::errstr;

my $sthnfo = $dbh->prepare( "SELECT * FROM SamInfo");
$sthnfo->execute();
my $qres = $sthnfo->fetchall_arrayref;
my %nfo;
if ($#$qres == 1) {
	for (@$qres) {
		my $Type = $$_[0];
		$nfo{"${Type}Ref"} = $_->[1];
		$nfo{"${Type}FQ1"} = $_->[2];
		$nfo{"${Type}FQ2"} = $_->[3];
		$nfo{"${Type}SAM"} = $_->[6];
	}
} else {die;}

#ddx \%nfo;
#   HumFQ1 => "s00_C.bs_1.fq.gz",
#   HumFQ2 => "s00_C.bs_2.fq.gz",
#   HumRef => "HomoGRCh38.fa",
#   HumSAM => "n3_grep.vircandi.sam.gz",
#   VirFQ1 => "n3_grep.1.fq.gz",
#   VirFQ2 => "n3_grep.2.fq.gz",
#   VirRef => "HBV.AJ507799.2.fa",
#   VirSAM => "n3_grep.vircandi.bshbv.bam",

getRefChrLen($nfo{'HumRef'}); getRefChrLen($nfo{'VirRef'});
#ddx \%ChrLen; ddx \%Genome;

my $sth = $dbh->prepare( "SELECT * FROM MergedSam
 WHERE HumChr IS NOT NULL AND VirChr IS NOT NULL AND HumCIAGR <> '*' AND VirCIAGR <> '*'
 ORDER BY HumChr,HumPos,VirChr,VirPos ASC;" );
$sth->execute();

# https://metacpan.org/pod/DBI#fetchall_arrayref
my $rows = []; # cache for batches of rows
while( my $row = ( shift(@$rows) || # get row from cache, or reload cache:
	shift(@{$rows=$sth->fetchall_arrayref(undef,10_000)||[]}) )
) {
	ddx $row;
=pod
# extractfuse.pl:28: [
#   "FCC1DVKACXX:7:1205:10704:40309#TTAGGCAT/2",
#   "AATAACNAAAAAAAAAAAAAAAAAAATAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAATAAATTAAACAAATAAAACAAAAAAAAAAA",
#   "bbbeeeBQ`ecgfhhiiiiiigea_\\TTGYbccaccOGJ]GJWOVaccccaa_ccccBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
#   163,
#   "chr11",
#   49559625,
#   "7S53M30S",
#   27,
#   163,
#   "gi|86261677|emb|AJ507799.2|",
#   96219,
#   "7S80M3S",
#   0,
#   "r,GA",
#   "r,GA",
# ]
=cut
	my $ret = Mgfq2HumVir($row,'Hum');
}

$dbh->rollback;
$dbh->disconnect;

__END__
../mergebam.pl n3_grep.vircandi.sam.gz n3_grep.vircandi.bshbv.bam n3_merged

../extractfuse.pl 200 n3_merged.sqlite n3
