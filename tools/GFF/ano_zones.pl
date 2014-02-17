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
#use Galaxy::ChromString;
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
    "dbi:SQLite:dbname=:memory:","","", 
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

warn "[!]Begin $annfile.\n";
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
#print "$id, $len\n";
}
#$ChrDat = ChrStrInit(\@IDs,\@Lens);
close A;

my %bitflag = (
		'CDS' => 1,
		'mRNA' => 2
	);

warn "[!]Begin $dbfile.\n";
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

sub check_mRNA_is_dup($$) {
	my ($mRNA_arr,$elm_arr)=@_;
	my ($mRNA_start,$mRNA_end,$mRNA_name)=@$mRNA_arr[2,3,6];
	for my $item (@$elm_arr) {
		if ( $$item[6] eq $mRNA_name ) {
			warn "ERR: CDS/UTR($$item[2]-$$item[3]) outsides mRNA:$mRNA_name($mRNA_start-$mRNA_end) !" if ( $mRNA_start > $$item[2] or $mRNA_end < $$item[3] );
			return 1;
		}
	}
	return 0;
}

warn "[!]Begin $infs.\n";
my @bamnames;
while (<$IN>) {
	if (/^#/) {
		if (/^#CHROM\tPOS\t/) {
			print O $_;
			#chomp;
			#@bamnames = split /\t/;
			#@bamnames = splice @bamnames, 9;
			#print O "# Bams: ",join(',',@bamnames),"\n";
#print "[@bamnames]\n";
		}
		next;
	}
	chomp;
	my ($chr, $pos, @data) = split /\t/; # undef, $refbase, $altbases, $QUAL, undef, $INFO, $FORMAT,
	if ($chr =~ /ref\|([^|]+)\|/) {
		$chr = $1;
	}
	$sthi->execute( $chr, $pos ) or die $sthi->errstr;
	my $qres = $sthi->fetchall_arrayref;
	if ($#$qres == -1) {
		print O join("\t", $chr, $pos, 'NA', @data),"\n";
		#print join("\t", $chr, $pos, 'NA'),"\n";
		next;
	}
	my (@CDS,@mRNA);
	for (@$qres) {
		if ($$_[1] =~ /mRNA/) {push @mRNA, join(',',@$_[6,2,3,4,5]);}
		 else {push @CDS, join(',',@$_[6,2,3,4,5]);}
	}
	if (@CDS > 0) {
		print O join("\t", $chr, $pos, 'CDS', @data, join(';',@CDS)),"\n";
		#print join("\t", $chr, $pos, 'CDS'),"\n";
		#ddx @CDS;
	} elsif (@mRNA > 0) {
		print O join("\t", $chr, $pos, 'mRNA', @data, join(';',@mRNA)),"\n";
		#print join("\t", $chr, $pos, 'mRNA'),"\n";
	} else { die; }
}
close $IN;

close O;

__END__
perl ano_zones.pl Dasnov3.db Dasnov3.ann bam2.bcgv.vcf.gz bam2.bcgv.vcf.ano

perl ano_zones.pl Dasnov3.db Dasnov3.ann bam2q20.snp.lst bam2q20.snp.ano
perl ano_zones.pl Dasnov3.db Dasnov3.ann bam2.bam.depth.gz bam2.bam.depth.ano

grep Total_Bases sai/*.fqstat | awk -F" " 'BEGIN {x=0} {x+=$2} END {print x}'
grep -P "\tNA\t" bam2q20.snp.ano | wc -l
grep -P "\tCDS\t" bam2q20.snp.ano | wc -l
grep -P "\tmRNA\t" bam2q20.snp.ano | wc -l

grep -P "\tNA\t" bam2.bam.depth.ano | wc -l
grep -P "\tNA\t" bam2.bam.depth.ano |  awk -F" " 'BEGIN {x=0} {x+=$4} END {print x}'
grep -P "\tCDS\t" bam2.bam.depth.ano | wc -l
grep -P "\tCDS\t" bam2.bam.depth.ano |  awk -F" " 'BEGIN {x=0} {x+=$4} END {print x}'
grep -P "\tmRNA\t" bam2.bam.depth.ano | wc -l
grep -P "\tmRNA\t" bam2.bam.depth.ano |  awk -F" " 'BEGIN {x=0} {x+=$4} END {print x}'
