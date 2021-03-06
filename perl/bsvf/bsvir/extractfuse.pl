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
$DEBUG=0.9;
our $minLen = 5;
our $low_complexity_cutoff_ratio=0.9; # 0.85
our $N_num_cutoff=2;

our (%Genome,%ChrLen);	# Let's share them.

my %attr = (
	RaiseError => 1,
	PrintError => 1,
	AutoCommit => 0
);
my $dbh = DBI->connect("dbi:SQLite:dbname=./$in","","",\%attr) or die $DBI::errstr;

unlink "./${out}.sqlite";
$dbh->sqlite_backup_to_file( "./$out.sqlite" );	# Well, we can just cp the file as well, before open it.
my $dboh = DBI->connect("dbi:SQLite:dbname=./${out}.sqlite","","",\%attr) or die $DBI::errstr;
$dboh->do("PRAGMA cache_size = -2000000");	# 2G http://www.sqlite.org/pragma.html#pragma_cache_size

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
 WHERE HumChr IS NOT NULL AND VirChr IS NOT NULL AND HumCIGAR <> '*' AND VirCIGAR <> '*'
 ORDER BY HumChr,HumPos,VirChr,VirPos ASC;" );
$sth->execute();

# https://metacpan.org/pod/DBI#fetchall_arrayref
my $rows = []; # cache for batches of rows
while( my $row = ( shift(@$rows) || # get row from cache, or reload cache:
	shift(@{$rows=$sth->fetchall_arrayref(undef,10_000)||[]}) )
) {
	#ddx $row;
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
	my $retHum = Mgfq2HumVir($row,'Hum');
	my $retVir = Mgfq2HumVir($row,'Vir');
	if ($#{$$retHum[0]}==1 and $#{$$retVir[0]}==1) {
		my ($PosBeginHum,$PosEndHum) = ( $$row[5]+$$retHum[0]->[0],$$row[5]+$$retHum[0]->[1] );
		my ($PosBeginVir,$PosEndVir) = ( $$row[10]+$$retVir[0]->[0],$$row[10]+$$retVir[0]->[1] );
		print join("\t",$$retHum[1],$$retVir[1],$$row[4],$PosBeginHum,$PosEndHum,$$row[9],$PosBeginVir,$PosEndVir),"\n";
	} elsif ($#{$$retHum[0]}==1 or $#{$$retVir[0]}==1) {
		if ( $DEBUG > 2 ) {
			ddx [$row,$retHum,$retVir];
			print "\nWARN\n\n";
		}
	}
}
$dbh->rollback;
$dbh->disconnect;

$dboh->commit;
$dboh->disconnect;
# 记录“3”的起止点。取1条，搜overlap，扩展start、end，再搜overlap，until差值恒定
# 一阶导数正负最大的就是两分界点。插入应该是酶反应，必然有重复序列。

__END__
../mergebam.pl n3_grep.vircandi.sam.gz n3_grep.vircandi.bshbv.bam n3_merged

../extractfuse.pl 200 n3_merged.sqlite n3
