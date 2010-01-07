#!/usr/bin/perl -w
use lib '/share/raid010/resequencing/resequencing/tmp/bin/annotation/glfsqlite';
use threads;
use strict;
use warnings;
use DBI;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.2.1;

our $opts='i:o:s:d:bvp';
our ($opt_i, $opt_o, $opt_s, $opt_v, $opt_b, $opt_d, $opt_p);

our $help=<<EOH;
\t-i Input SNP file (Human.Q20)
\t-p Input SNP file is of a population (*.add_ref)
\t  [Chromosome name will be striped.]
\t-s Specie name (human)
\t  [Specie name MUST be ths SAME throughout !]
\t-d Genome SQLite data file (_dbGFF.sqlite)
\t-o Output SQLite data file (_result.sqlite)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_i='Human.Q20' if ! defined $opt_i;
$opt_o='_result.sqlite' if ! $opt_o;
$opt_d='_dbGFF.sqlite' if ! $opt_d;
$opt_s='human' if ! $opt_s;

print STDERR "From [$opt_i][$opt_d] to [$opt_o], Specie:[$opt_s]\n";
if ($opt_p) {print STDERR "Population mode ON.\n";}
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];

my %attr = (
    RaiseError => 0,
    PrintError => 1,
    AutoCommit => 0
);

### /dev/shm
my $shm_real='/dev/shm/sqlite_mirror.'.$$;
my $shm_real_in='/dev/shm/sqlite_in.'.$$;
unlink ($shm_real,$shm_real_in);	# Well, what if the computer rebooted and you are so lucky ...
system 'cp','-pf',$opt_d,$shm_real_in;
###

my $dbh = DBI->connect('dbi:SQLite:dbname='.$shm_real_in,'','',\%attr) or die $DBI::errstr;
our $rdbh = DBI->connect('dbi:SQLite:dbname='.$shm_real,'','',\%attr) or die $DBI::errstr;

my $sql=q/
CREATE TABLE IF NOT EXISTS res{---}
(
   chrid TEXT,
   position INTEGER,
   refseq TEXT,
   snpseq TEXT,
   snp_q INTEGER,
   rna_chg TEXT,
   aa_chg TEXT,
   chged INTEGER,
   primary_inf TEXT,
   start INTEGER,
   end INTEGER,
   strand TEXT,
   frame TEXT,
   name TEXT,
   groups TEXT
);
/;
for (split /;/,$sql) {
	next if /^\s*$/;
	s/{---}/$opt_s/g;
	$rdbh->do($_) or die $rdbh->errstr;
}
$rdbh->commit;
my $sthi = $dbh->prepare( "SELECT * FROM gff$opt_s WHERE chrid=? AND ? BETWEEN start AND end" );
my $rsthi = $rdbh->prepare( "INSERT INTO res$opt_s ( chrid,position,refseq,snpseq,snp_q,primary_inf,start,end,strand,frame,name,groups ) VALUES ( ?,?,?,?,?,?,?,?,?,?,?,? )" );

 my %IUB = ( A => [qw(A)],
	     C => [qw(C)],
	     G => [qw(G)],
	     T => [qw(T)],
	     U => [qw(U)],
	     M => [qw(A C)],
	     R => [qw(A G)],
	     W => [qw(A T)],
	     S => [qw(C G)],
	     Y => [qw(C T)],
	     K => [qw(G T)],
	     V => [qw(A C G)],
	     H => [qw(A C T)],
	     D => [qw(A G T)],
	     B => [qw(C G T)],
	     X => [qw(G A T C)],
	     N => [qw(G A T C)]
	     );

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

sub do_snp($$$) {
    my $file=$_[0];
    my $dbhs=$_[1];
    my $rdbhs=$_[2];
    die "FATAL: Need a SCALAR reference !\n" if ref($file) ne 'SCALAR';
    my (@dat,@snpseqd,$qres,$snpseqds,$aa_chg);
    open FILE,'<',$$file or die "Error: $!\n";
    while (<FILE>) {
	next if /^\s+$/ or /^;+/;
	chomp;
	my ($seqname, $position, $refseq, $snpseq, $snp_q);
	@snpseqd=();
	$snpseqds='';
	if ($opt_p) {
		my ($snp_pop,@snpop,%bases);
		($seqname, $position, $refseq, $snp_pop) = split /\t/;
		$snp_q=100;	# dummy
		@snpop=split / /,$snp_pop;
		for (@snpop) {
			next if $_ eq '-';
			$snpseq=$IUB{$_};
			for (@$snpseq) {
				next if $_ eq $refseq;
				++$bases{$_};
			}
		}
		@snpseqd=keys %bases;
		$snpseqds=join '',@snpseqd;
	} else {
		($seqname, $position, $refseq, $snpseq, $snp_q) = split /\t/;
		$snpseq=$IUB{$snpseq};
		for (@$snpseq) {
		   next if $_ eq $refseq;
		   push @snpseqd,$_;
		   $snpseqds .= $_;
		}
		warn "Not 1 SNP exists as $seqname\t$position\t$refseq\t@${snpseq}\t${snp_q} !\n" if $#snpseqd!=0;
	}
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
	$$dbhs->execute( $seqname, $position ) or die $$dbhs->errstr;
	$qres = $$dbhs->fetchall_arrayref;
	if ($#$qres == -1) {
	   next;
	}
	if ($#$qres == 0) {
	   my $res=$$qres[0];
	   for (@snpseqd) {
		$$rdbhs->execute( $seqname,$position,$refseq,$_,$snp_q,
		   @$res[1,2,3,4,5,6,7] ) or die $$rdbhs->errstr;
		print "$seqname,$position,$refseq,$_,$snp_q,@$res[1,2,3,4,5,6]\n" if $opt_v;
	   }
	   next;
	}
	my (@Elements,@mRNA)=();
	for (@$qres) {
		if ($$_[1] =~ /mRNA/) {push @mRNA,$_;}
		 else {push @Elements,$_;}
	}
	for my $mRNA_arr (@mRNA) {
		next if ! defined $mRNA_arr;
		if ( check_mRNA_is_dup($mRNA_arr,\@Elements) ) { $mRNA_arr = undef; next;}
	}
	# now, no duped mRNA.
	for my $res (@Elements,@mRNA) {
	   next if ! defined $res;
	   for (@snpseqd) {
		$$rdbhs->execute( $seqname,$position,$refseq,$_,$snp_q,
		   @$res[1,2,3,4,5,6,7] )
		   or die $$rdbhs->errstr;
		print "$seqname,$position,$refseq,$_,$snp_q,@$res[1,2,3,4,5,6,7]\n" if $opt_v;
	   }
	}
    }
    $$dbhs->finish;
    $$rdbhs->finish;
    close FILE;
    return 1;
}

my $name=$opt_i;
do_snp(\$name,\$sthi,\$rsthi);
$rdbh->commit;

$dbh->rollback or warn $dbh->errstr;
$dbh->disconnect;

$sql=q/
CREATE INDEX IF NOT EXISTS ncs{---} ON res{---}(name,chrid,position);
CREATE INDEX IF NOT EXISTS pi{---} ON res{---}(primary_inf);
CREATE INDEX IF NOT EXISTS pi{---} ON res{---}(chged);
/;
for (split /;/,$sql) {
	next if /^\s*$/;
	s/{---}/$opt_s/g;
	$rdbh->do($_) or warn $rdbh->errstr;
}
$rdbh->commit;
$rdbh->disconnect;

my $read_time = [gettimeofday];
my $thr1 = async { system 'cp','-pf',$shm_real,$opt_o; };
my $thr2 = async {
	system 'bzip2','-9k',$shm_real;
	system 'mv','-f',$shm_real.'.bz2',$opt_o.'.bz2';
};
$thr1->join();
my $copy_time = [gettimeofday];
$thr2->join();
unlink $shm_real.'.bz2';
unlink ($shm_real,$shm_real_in);
###

my $stop_time = [gettimeofday];

$|=1;
print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time ),
	" second(s).\n   Parseing file used\t",tv_interval( $start_time, $read_time ),
	" second(s).\n   Moving SQLite file used\t",tv_interval( $read_time, $copy_time )," second(s).\n";

print STDERR "\033[32;1m Please use [$opt_s] as Specie name in later steps.\033[0;0m\n";
