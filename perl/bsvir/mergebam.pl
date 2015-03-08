#!/usr/bin/env perl
use strict;
use warnings;
use DBI qw(:sql_types);
use Data::Dump qw(ddx);

die "Usage: $0 <hum_sorted_bam> <vir_sorted_bam> <out_prefix>\n" if @ARGV < 3;
my ($inhum,$invir,$out)=@ARGV;

unlink "./${out}.sqlite";
my $dbh = DBI->connect("dbi:SQLite:dbname=./${out}.sqlite","","",{
	RaiseError => 0,
	PrintError => 1,
	AutoCommit => 0
}) or die $DBI::errstr;
$dbh->do("PRAGMA cache_size = -2000000");	# 2G http://www.sqlite.org/pragma.html#pragma_cache_size
print '[!]SQLite version: ',$dbh->{sqlite_version}," (3.7.10 or above is better).\n";

my $sql=q/
CREATE TABLE MergedSam
(  FQID12 TEXT,
   FQSeq TEXT,
   FQQual TEXT,
   HumFlag INTEGER,
   HumChr TEXT,
   HumPos INTEGER,
   HumCIAGR TEXT,
   HumMapQ INTEGER,
   VirFlag INTEGER,
   VirChr TEXT,
   VirPos INTEGER,
   VirCIAGR TEXT,
   VirMapQ INTEGER,
   YDYCHum TEXT,
   YDYCVir TEXT
);
CREATE INDEX nFQID12 ON MergedSam(FQID12);
/;
for (split /;\s*/,$sql) {
	next if /^\s*$/;
	$dbh->do($_) or die $dbh->errstr,"[$_]";
}
$dbh->commit;

my $inStr0 = "INSERT INTO MergedSam (FQID12,FQSeq,FQQual,{-}Flag,{-}Chr,{-}Pos,{-}CIAGR,{-}MapQ,YDYC{-}) VALUES (?,?,?,?,?,?,?,?,?)";
my ($inStr1,$inStr2) = ($inStr0,$inStr0); $inStr1 =~ s/{-}/Hum/g; $inStr2 =~ s/{-}/Vir/g;
my $sthins1 = $dbh->prepare($inStr1); my $sthins2 = $dbh->prepare($inStr2);

my $upStr0 = "UPDATE MergedSam SET {-}Flag=?,{-}Chr=?,{-}Pos=?,{-}CIAGR=?,{-}MapQ=?,YDYC{-}=? WHERE FQID12=?";
my ($upStr1,$upStr2) = ($upStr0,$upStr0); $upStr1 =~ s/{-}/Hum/g; $upStr2 =~ s/{-}/Vir/g;
my $sthupd1 = $dbh->prepare($upStr1); my $sthupd2 = $dbh->prepare($upStr2);

my $sthid = $dbh->prepare( "SELECT * FROM MergedSam WHERE FQID12=?" );

sub Sam2FQ($) {
	my $rd = $_[0];
	my ($rd12);
	if ($$rd[1] & 0x40) {
		$rd12 = 1;
	} elsif ($$rd[1] & 0x80) {
		$rd12 = 2;
	}
	die "CIAGR:$$rd[5]" if $$rd[5] =~ /H/;
	my $seq = $$rd[9];
	my $qual = $$rd[10];
	if ($$rd[1] & 0x10) {
		$seq =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
		$seq = reverse $seq;
		$qual = reverse $qual;
	}
	my $id = join('',$$rd[0],'/',$rd12);
	return [$rd12,$id,$seq,$qual];
}
sub InsertSam($$) {
	my ($Type,$in)=@_;
	open( IN,"-|","samtools view -F768 $in") or die "Error opening $in: $!\n";
	while (my $line = <IN>) {
		#my ($id, $flag, $ref, $pos, $mapq, $CIAGR, $mref, $mpos, $isize, $seq, $qual, @OPT) = split /\t/,$line;
		#print "$id, $flag, $ref, $pos, $mapq, $CIAGR, $mref, $mpos, $isize\n";
		my @Dat1 = split /\t/,$line;
		my ($id, $flag, $ref, $pos, $mapq, $CIAGR, $mref, $mpos, $isize, $seq, $qual, @OPT) = @Dat1;
		my $ret1 = Sam2FQ(\@Dat1);
		(undef,$id,$seq,$qual) = @$ret1;
		my $YDYC;
		{
			my ($YD,$YC)=('','');
			for (@OPT) {
				if (/YD:Z:(\w+)/) {
					$YD = $1;
				}
				$YC=$1 if /YC:Z:(\w+)/
			}
			$YDYC = "$YD,$YC";
		}
		$sthid->execute($id);
		my $qres = $sthid->fetchall_arrayref;
		if ($#$qres == -1) {
			$sthins1->execute($id,$seq,$qual,$flag,$ref,$pos,$CIAGR,$mapq,$YDYC) if $Type eq 'Hum';
			$sthins2->execute($id,$seq,$qual,$flag,$ref,$pos,$CIAGR,$mapq,$YDYC) if $Type eq 'Vir';
		} elsif ($#$qres == 0) {
			$sthupd1->execute($flag,$ref,$pos,$CIAGR,$mapq,$YDYC,$id) if $Type eq 'Hum';
			$sthupd2->execute($flag,$ref,$pos,$CIAGR,$mapq,$YDYC,$id) if $Type eq 'Vir';
		} else {die;}
	}
}

InsertSam('Hum',$inhum);
$dbh->commit;
InsertSam('Vir',$invir);
$dbh->commit;





$dbh->commit;
$dbh->disconnect;

__END__
../mergebam.pl n3_grep.vircandi.sam.gz n3_grep.vircandi.bshbv.bam n3_merged
sqlite3 n3_merged.sqlite .dump >n3_merged.sqlite.dump
