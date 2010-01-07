#!/usr/bin/perl -w
use lib '/share/raid010/resequencing/resequencing/tmp/bin/annotation/glfsqlite';
use threads;
use strict;
use warnings;
use DBI;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.2.2;

our $opts='i:o:s:d:bvf';
our ($opt_i, $opt_o, $opt_s, $opt_v, $opt_b, $opt_d, $opt_f);

our $help=<<EOH;
\t-i Input Genome sequence file (human.fa)
\t  [Chromosome name will be striped.]
\t-s Specie of the GFF3 (human)
\t  [Specie name MUST be ths SAME throughout !]
\t-d Genome SQLite data file (_dbGFF.sqlite)
\t-f fix for wrong GFF frame
\t-o Output SQLite data file (_result.sqlite)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_i='human.fa' if ! defined $opt_i;
$opt_o='_result.sqlite' if ! $opt_o;
$opt_d='_dbGFF.sqlite' if ! $opt_d;
$opt_s='human' if ! $opt_s;

print STDERR "Going to alter GFF frame for BGI version GFF.\n" if $opt_f;
print STDERR "From [$opt_i][$opt_d] to update [$opt_o], Specie:[$opt_s]\n";
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];

my $shm_real='/dev/shm/sqlite_mirror.'.$$;
system 'cp','-pf',$opt_o,$shm_real;
my $shm_real_in='/dev/shm/sqlite_in.'.$$;
system 'cp','-pf',$opt_d,$shm_real_in;

my %attr = (
    RaiseError => 0,
    PrintError => 1,
    AutoCommit => 0
);
my $dbh = DBI->connect('dbi:SQLite:dbname='.$shm_real_in,'','',\%attr) or die $DBI::errstr;
our $rdbh = DBI->connect('dbi:SQLite:dbname='.$shm_real,'','',\%attr) or die $DBI::errstr;
my $sql=q/
CREATE INDEX IF NOT EXISTS ncs{---} ON res{---}(name,chrid,position);
CREATE INDEX IF NOT EXISTS pi{---} ON res{---}(primary_inf);
CREATE INDEX IF NOT EXISTS pi{---} ON res{---}(chged);
/;
for (split /;/,$sql) {
	next if /^\s*$/;
	s/{---}/$opt_s/g;
	$rdbh->do($_) or die $rdbh->errstr;
}
###################
our ($sdh,$sdh2,$srdhi,$srdho);	# our so that to be used in subs
sub prepareSQL($) {
	my $specname=$_[0];
	my $sql=q/
SELECT * FROM gff/.$specname.q/ WHERE chrid=? AND name=? AND primary_inf LIKE '%mRNA'
 AND ? BETWEEN start AND end
/;
	$sdh = $dbh->prepare($sql) or warn $dbh->errstr;
	$sql=q/
SELECT * FROM gff/.$specname.q/ WHERE chrid=:1 AND name=:2 AND primary_inf LIKE '%CDS'
 AND start BETWEEN :3 AND :4
 AND end BETWEEN :3 AND :4
 ORDER BY start ASC
/;
	$sdh2 = $dbh->prepare($sql) or warn $dbh->errstr;
#########
	$sql=q/
SELECT chrid,position,refseq,snpseq,primary_inf,name FROM res/.$specname.q/ WHERE aa_chg IS NULL AND primary_inf LIKE '%CDS'
/;
	$srdhi = $rdbh->prepare($sql) or warn $rdbh->errstr;
	$srdhi->execute;
	$sql=q/
UPDATE res/.$specname.q/ SET rna_chg=?,aa_chg=?,chged=? WHERE name=? AND chrid=? AND position=?
/;
	$srdho = $rdbh->prepare($sql) or warn $rdbh->errstr;
}
###################
my %gen_code = ( "TTT" => "Phe", "TTC" => "Phe", "TTA" => "Leu", "TTG" => "Leu", "CTT" => "Leu", "CTC" => "Leu", "CTA" => "Leu", "CTG" => "Leu", "ATT" => "Ile", "ATC" => "Ile", "ATA" => "Ile", "ATG" => "Met", "GTT" => "Val", "GTC" => "Val", "GTA" => "Val", "GTG" => "Val", "TCT" => "Ser", "TCC" => "Ser", "TCA" => "Ser", "TCG" => "Ser", "CCT" => "Pro", "CCC" => "Pro", "CCA" => "Pro", "CCG" => "Pro", "ACT" => "Thr", "ACC" => "Thr", "ACA" => "Thr", "ACG" => "Thr", "GCT" => "Ala", "GCC" => "Ala", "GCA" => "Ala", "GCG" => "Ala",  "TAT" => "Tyr", "TAC" => "Tyr", "TAA" => "STOP", "TAG" => "STOP", "CAT" => "His", "CAC" => "His", "CAA" => "Gln", "CAG" => "Gln", "AAT" => "Asn", "AAC" => "Asn", "AAA" => "Lys", "AAG" => "Lys", "GAT" => "Asp", "GAC" => "Asp", "GAA" => "Glu", "GAG" => "Glu",  "TGT" => "Cys", "TGC" => "Cys", "TGA" => "STOP", "TGG" => "Trp", "CGT" => "Arg", "CGC" => "Arg", "CGA" => "Arg", "CGG" => "Arg", "AGT" => "Ser", "AGC" => "Ser", "AGA" => "Arg", "AGG" => "Arg", "GGT" => "Gly", "GGC" => "Gly", "GGA" => "Gly", "GGG" => "Gly" );
###################
sub q_3CDS($$$) {
	my ($seqname,$position,$name)=@_;
	my ($mRNAstart,$mRNAend,$mRNAstrand,@CDS,$CDSpointer);
	$sdh->execute( ($seqname,$name,$position) );
	my $qres = $sdh->fetchall_arrayref;
	if ($#$qres != 0) {
		if ($#$qres == -1) {
			warn "No info. for $name\t$seqname\t$position !\n";
			return [[],'NA'];
		} else {warn "$#$qres more hit(s) for $name\$seqname\t$position !\n";}
	}
	$mRNAstart=$$qres[0][2];
	$mRNAend=$$qres[0][3];
	$mRNAstrand=$$qres[0][4];
	@CDS=();
	$sdh2->execute( ($seqname,$name,$mRNAstart,$mRNAend) );
	$qres = $sdh2->fetchall_arrayref;
	if ($#$qres == -1) {
		warn "No info. for $seqname\t$mRNAstart-$mRNAend !\n";
	}
	for (@$qres) {
		warn "$$_[4] NOT the same as $mRNAstrand for  $$_[6]($$_[0] $$_[1] $$_[2] $$_[3])\n"
		   if $$_[4] ne $mRNAstrand;
		if ($opt_f) {
			if ($$_[5]==1) {$$_[5]=2}
			 elsif ($$_[5]==2) {$$_[5]=1}
		}
		push @CDS,[$$_[2],$$_[3],$$_[5]];	# start, end, frame
		}

	$CDSpointer=-1;	# well, pop will be faster, but, pointer is easy to write.
	for (@CDS) {
		++$CDSpointer;
		if ($position>$$_[1]) {
			next;
		} elsif ($position>=$$_[0]) {
			last;
		}
	}
	if ($CDSpointer>0) {
		@CDS=@CDS[$CDSpointer-1,$CDSpointer,$CDSpointer+1];
	} else {@CDS=(undef,@CDS[$CDSpointer,$CDSpointer+1]);}
	return [\@CDS,$mRNAstrand];
}
###################
sub q_subseq($$$$$) {
	my ($genome_hash,$seqname,$c_CDS_arr,$strand,$position)=@_;	# all starting 1
	my ($forward_strain,$sense_strain,$pos_rel);
	$forward_strain=substr $$genome_hash{$seqname},$$c_CDS_arr[0]-1,$$c_CDS_arr[1]-$$c_CDS_arr[0]+1;
	if ($strand eq '+') {
		$pos_rel=$position-$$c_CDS_arr[0];	# +1 for starting 1, so no chg for starting 0
		$sense_strain=$forward_strain;
	} else {
		$pos_rel=$$c_CDS_arr[1]-$position;	# +1 for starting 1, so no chg for starting 0
		$forward_strain =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;	# http://doc.bioperl.org/releases/bioperl-1.4/Bio/Tools/SeqPattern.html#CODE6
		$sense_strain=reverse $forward_strain;
	}
	return [\$sense_strain,$pos_rel];
}
###################
sub q_aa_chg($$$$$) {
	my ($sense_seq_ref,$pos_rel,$ref_base,$snp_base,$strand)=@_;	# $pos_rel starting 0
	if ($strand eq '-') {
		$ref_base =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
		$snp_base =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
	}
	my $gom_base=substr $$sense_seq_ref,$pos_rel,1;
	$ref_base = uc $ref_base; $gom_base = uc $gom_base;
	warn "WARING: ref_base($ref_base)[$strand] differ to that($gom_base) in Genome [$$sense_seq_ref]!\n" if $ref_base ne $gom_base;	# well, let's smoking...
	my $posIN3=$pos_rel % 3;
	my $codon_start = $pos_rel-$posIN3;
	my $ref_seq=substr $$sense_seq_ref,$codon_start,3;
	my $snp_seq=substr($ref_seq,0,$posIN3).$snp_base.substr($ref_seq,$posIN3+1,2-$posIN3);
	$ref_seq = uc $ref_seq; $snp_seq = uc $snp_seq;
	my ($aa_ref,$aa_snp)=($gen_code{$ref_seq},$gen_code{$snp_seq});
	unless ( defined $aa_ref or defined $aa_snp ) {
		warn "ERROR: [$ref_seq] or [$snp_seq] not translatable !\n";
		return ["$ref_seq->$snp_seq",'Unknown',-1];
	}
	if ($aa_ref eq $aa_snp) {
		return ["$ref_seq->$snp_seq","$aa_ref->$aa_snp",0];
	} else {return ["$ref_seq->$snp_seq","$aa_ref->$aa_snp",1];}
}
###################
sub q_aa($$$$$$$) {
	my ($genome_hash,$seqname,$CDS_arr,$strand,$position,$ref_base,$snp_base)=@_;
	my (@centerCDS,$aa_chg,$sense_seq_ref,$pos_rel,$comp_seq_ref,@compCDS);
	if ($strand eq '+') {	# @$CDS_arr : start, end, frame (starting 1)
		@centerCDS=( $CDS_arr->[1]->[0]+$CDS_arr->[1]->[2],$CDS_arr->[1]->[1] );
		@compCDS = ( $CDS_arr->[2]->[0],$CDS_arr->[2]->[0]+1 ) if defined($CDS_arr->[2]) and $CDS_arr->[2]->[0]+1 <= $CDS_arr->[2]->[1];
	} elsif ($strand eq '-') {
		@centerCDS=( $CDS_arr->[1]->[0],$CDS_arr->[1]->[1]-$CDS_arr->[1]->[2] );
		@compCDS = ( $CDS_arr->[0]->[1]-1,$CDS_arr->[0]->[1] ) if defined($CDS_arr->[0]) and $CDS_arr->[0]->[1]-1 >= $CDS_arr->[0]->[0];
	} else { return []; }	# $strand eq 'NA'
	if ($position >= $centerCDS[0] and $position <= $centerCDS[1]) {
		($sense_seq_ref,$pos_rel)=@{q_subseq($genome_hash,$seqname,\@centerCDS,$strand,$position)};
		if ( $#compCDS == 1 ) {
			($comp_seq_ref,undef)=@{q_subseq($genome_hash,$seqname,\@compCDS,$strand,$position)};
			$$comp_seq_ref = $$sense_seq_ref.$$comp_seq_ref;
		} else {$$comp_seq_ref = $$sense_seq_ref;}
	} elsif ( $position < $centerCDS[0] and $strand eq '+' ) {
		($sense_seq_ref,$pos_rel)=@{q_subseq($genome_hash,$seqname,$CDS_arr->[1],$strand,$position)};
		@compCDS=( $CDS_arr->[0]->[1] - 2+$CDS_arr->[1]->[2],$CDS_arr->[0]->[1] );
		($comp_seq_ref,undef)=@{q_subseq($genome_hash,$seqname,\@compCDS,$strand,$position)};
		$$comp_seq_ref .= $$sense_seq_ref;
		$pos_rel += 3-$CDS_arr->[1]->[2];
	} elsif ( $position > $centerCDS[0] and $strand eq '-' ) {
		($sense_seq_ref,$pos_rel)=@{q_subseq($genome_hash,$seqname,$$CDS_arr[1],$strand,$position)};
		@compCDS=( $CDS_arr->[2]->[0],$CDS_arr->[2]->[0] + 2-$CDS_arr->[1]->[1] );
		($comp_seq_ref,undef)=@{q_subseq($genome_hash,$seqname,\@compCDS,$strand,$position)};
		$$comp_seq_ref .= $$sense_seq_ref;
		$pos_rel += 3-$CDS_arr->[1]->[2];
	} else {
		warn "ERROR in q_aa() !\n$position,[$CDS_arr->[1]->[0],$CDS_arr->[1]->[1]] $strand ";
		return [];	# well, a bit smoke
	}
	my @aa_chg=@{q_aa_chg($comp_seq_ref,$pos_rel,$ref_base,$snp_base,$strand)};
	return \@aa_chg;
}
###################
sub load_genome($) {
	open IN,'<',$_[0] or die "Error: $!\n";
	my %genome;
	while (<IN>) {
		s/^>//;
		my $title = $_;
		my $seqname = $1 if($title =~ /^(\S+)/);
		$seqname =~ s/^chr
			(?>
				((?<=^chr)o)?
				((?<=^chro)m)?
				((?<=^chrom)o)?
				((?<=^chromo)s)?
				((?<=^chromos)o)?
				((?<=^chromoso)m)?
				((?<=^chromosom)e)?
			)//xi;
warn "loading > $seqname ...\n";
		$/=">";
		$genome{$seqname}=<IN>;
		chomp $genome{$seqname};
		$genome{$seqname}=~s/\s//g;
		$/="\n";
	}
	close IN;
	return \%genome;
}
#main()
my $genome_hash;
$genome_hash=load_genome($opt_i);
prepareSQL($opt_s);	# No time to rewrite now ...
while (my $ary_ref = $srdhi->fetchrow_arrayref) {
	my ($chrid,$position,$ref_base,$snp_base,$primary_inf,$name)=@$ary_ref;
	my $rv=q_3CDS($chrid,$position,$name);	# [\@CDS,$mRNAstrand]
	my $chg_ref=q_aa( $genome_hash,$chrid,$$rv[0],$$rv[1],$position,$ref_base,$snp_base );
	next if $#$chg_ref == -1;
	my ($rna_chg,$aa_chg,$chged)=@$chg_ref;
print "$name\t$chrid\t$position\t$ref_base -> $snp_base\t$rna_chg\t$aa_chg\t$chged\n" if $opt_v;
	$srdho->execute($rna_chg,$aa_chg,$chged,$name,$chrid,$position) or warn $rdbh->errstr;
}
$rdbh->commit or warn $rdbh->errstr;
$rdbh->disconnect;
$dbh->rollback or warn $dbh->errstr;
$dbh->disconnect;

my $read_time = [gettimeofday];
my $thr1 = async { system 'cp','-pf',$shm_real,$opt_o; };
my $thr2 = async {
	system 'bzip2','-9k',$shm_real;
	system 'mv','-f',$shm_real.'.bz2',$opt_o.'.bz2';
};
$thr1->join();
my $copy_time = [gettimeofday];
$thr2->join();
unlink ($shm_real,$shm_real_in);
unlink $shm_real.'.bz2';

my $stop_time = [gettimeofday];

$|=1;
print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time ),
	" second(s).\n   Parseing file used\t",tv_interval( $start_time, $read_time ),
	" second(s).\n   Moving SQLite file used\t",tv_interval( $read_time, $copy_time )," second(s).\n";

print STDERR "\033[32;1m Please use [$opt_s] as Specie name in later steps.\033[0;0m\n";
__END__

./update_aa.pl -b -o _result1.sqlite > _result1.txt 2>_result1.err &
