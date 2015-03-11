#!/usr/bin/env perl
use strict;
use warnings;
use DBI qw(:sql_types);
use File::Basename;
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
   HumCIGAR TEXT,
   HumMapQ INTEGER,
   VirFlag INTEGER,
   VirChr TEXT,
   VirPos INTEGER,
   VirCIGAR TEXT,
   VirMapQ INTEGER,
   YDYCHum TEXT,
   YDYCVir TEXT
);
CREATE INDEX nFQID12 ON MergedSam(FQID12);
CREATE TABLE SamInfo (
   Type TEXT,
   Ref TEXT,
   FQ1 TEXT,
   FQ2 TEXT,
   ReadLen1 INTEGER,
   ReadLen2 INTEGER,
   Sam TEXT
);
/;
for (split /;\s*/,$sql) {
	next if /^\s*$/;
	$dbh->do($_) or die $dbh->errstr,"[$_]";
}
$dbh->commit;

my $inStr0 = "INSERT INTO MergedSam (FQID12,FQSeq,FQQual,{-}Flag,{-}Chr,{-}Pos,{-}CIGAR,{-}MapQ,YDYC{-}) VALUES (?,?,?,?,?,?,?,?,?)";
my ($inStr1,$inStr2) = ($inStr0,$inStr0); $inStr1 =~ s/{-}/Hum/g; $inStr2 =~ s/{-}/Vir/g;
my $sthins1 = $dbh->prepare($inStr1); my $sthins2 = $dbh->prepare($inStr2);

my $upStr0 = "UPDATE MergedSam SET {-}Flag=?,{-}Chr=?,{-}Pos=?,{-}CIGAR=?,{-}MapQ=?,YDYC{-}=? WHERE FQID12=?";
my ($upStr1,$upStr2) = ($upStr0,$upStr0); $upStr1 =~ s/{-}/Hum/g; $upStr2 =~ s/{-}/Vir/g;
my $sthupd1 = $dbh->prepare($upStr1); my $sthupd2 = $dbh->prepare($upStr2);

my $sthid = $dbh->prepare( "SELECT * FROM MergedSam WHERE FQID12=?" );

sub openfile($) {
	my ($filename)=@_;
	my $infile;
	if ($filename=~/.xz$/) {
		open( $infile,"-|","xz -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.gz$/) {
		open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.bz2$/) {
		open( $infile,"-|","bzip2 -dc $filename") or die "Error opening $filename: $!\n";
	} else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
	return $infile;
}
sub getsamChrLen($$) {
	my ($Type,$in) = @_;
	my $inpath = dirname($in);
	my ($Ref,$fq1,$fq2,$readlen1,$readlen2);
	open( IN,"-|","samtools view -H $in") or die "Error(m1) opening $in: $!\n";
	while (<IN>) {
		chomp;
		my ($t,$id,$ln)=split /\t/;
		if ($t eq '@SQ') {
			$id = (split /\:/,$id)[1];
			$ln = (split /\:/,$ln)[1];
			#$ChrLen{$id} = $ln;
		} elsif ($t eq '@PG') {
# @PG	ID:BSMAP	VN:2.87	CL:"bsmap -u -d HomoGRCh38.fa -a F12HPCCCSZ0010_Upload/s00_C.bs_1.fq.gz -b F12HPCCCSZ0010_Upload/s00_C.bs_2.fq.gz -z 64 -p 12 -v 10 -q 2 -o s00_C.bs.bam"
# @PG	ID:bwa-meth	PN:bwa-meth	VN:0.10	CL:"./bwameth.py -t 24 -p s00_C.bshum --read-group s00_C --reference HomoGRCh38.fa s00_C.bs_1.fq.gz s00_C.bs_2.fq.gz"
			if ($id eq 'ID:BSMAP') {
				$Ref = $1 if /\-d\s+([.\w]+)\b/;
				$fq1 = $1 if /\-a\s+([^\s"]+)(\s|"?$)/;
				$fq2 = $1 if /\-b\s+([^\s"]+)(\s|"?$)/;
			} elsif ($id eq 'ID:bwa-meth') {
				$Ref = $1 if /\--reference\s+([.\w]+)\b/;
				if (/CL:.+\s([^-\s"]+)\s+([^-\s"]+)(\s|"?$)/) {
					$fq1 = $1;
					$fq2 = $2;
				}
			} else {die;}

		}
	}
	if ($inpath ne '.') {
		($Ref,$fq1,$fq2) = map { "$inpath/$_" } ($Ref,$fq1,$fq2);
	}
#warn "$Ref,$fq1,$fq2,$inpath";
	my $FQ1 = openfile($fq1);
	<$FQ1>;chomp($_=<$FQ1>);
	$readlen1 = length $_;
	close $FQ1;
	my $FQ2 = openfile($fq2);
	<$FQ2>;chomp($_=<$FQ2>);
	$readlen2 = length $_;
	close $FQ2;
	my $sthinsnfo = $dbh->prepare("INSERT INTO SamInfo (Type,Ref,FQ1,FQ2,ReadLen1,ReadLen2,Sam) VALUES (?,?,?,?,?,?,?)");
	$sthinsnfo->execute($Type,$Ref,$fq1,$fq2,$readlen1,$readlen2,$in);
	warn "[!]Reading [$in]: ($Type,$Ref,$fq1,$fq2,$readlen1,$readlen2)\n";
}

sub Sam2FQ($) {
	my $rd = $_[0];
	my ($rd12);
	if ($$rd[1] & 0x40) {
		$rd12 = 1;
	} elsif ($$rd[1] & 0x80) {
		$rd12 = 2;
	}
	die "CIGAR:$$rd[5]" if $$rd[5] =~ /H/;
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
	getsamChrLen($Type,$in);
	open( IN,"-|","samtools view -F768 $in") or die "Error opening $in: $!\n";
	while (my $line = <IN>) {
		#my ($id, $flag, $ref, $pos, $mapq, $CIGAR, $mref, $mpos, $isize, $seq, $qual, @OPT) = split /\t/,$line;
		#print "$id, $flag, $ref, $pos, $mapq, $CIGAR, $mref, $mpos, $isize\n";
		my @Dat1 = split /\t/,$line;
		my ($id, $flag, $ref, $pos, $mapq, $CIGAR, $mref, $mpos, $isize, $seq, $qual, @OPT) = @Dat1;
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
			$sthins1->execute($id,$seq,$qual,$flag,$ref,$pos,$CIGAR,$mapq,$YDYC) if $Type eq 'Hum';
			$sthins2->execute($id,$seq,$qual,$flag,$ref,$pos,$CIGAR,$mapq,$YDYC) if $Type eq 'Vir';
		} elsif ($#$qres == 0) {
			$sthupd1->execute($flag,$ref,$pos,$CIGAR,$mapq,$YDYC,$id) if $Type eq 'Hum';
			$sthupd2->execute($flag,$ref,$pos,$CIGAR,$mapq,$YDYC,$id) if $Type eq 'Vir';
		} else {die;}
	}
}

InsertSam('Hum',$inhum);
$dbh->commit;
InsertSam('Vir',$invir);
$dbh->commit;



$sql=q/
CREATE INDEX nSort ON MergedSam(HumChr,HumPos,VirChr,VirPos);
/;
for (split /;\s*/,$sql) {
	next if /^\s*$/;
	$dbh->do($_) or die $dbh->errstr,"[$_]";
}
$dbh->commit;
$dbh->disconnect;

__END__
../mergebam.pl n3_grep.vircandi.sam.gz n3_grep.vircandi.bshbv.bam n3_merged
sqlite3 n3_merged.sqlite .dump >n3_merged.sqlite.dump


sqlite> EXPLAIN QUERY PLAN SELECT * FROM MergedSam WHERE HumChr IS NOT NULL AND VirChr IS NOT NULL AND HumCIGAR <> '*' AND VirCIGAR <> '*' ORDER BY HumChr,HumPos,VirChr,VirPos ASC;
0|0|0|SCAN TABLE MergedSam USING INDEX nSort (~62500 rows)
