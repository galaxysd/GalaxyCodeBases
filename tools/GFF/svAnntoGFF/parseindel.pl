#!/bin/env perl
use lib '/share/raid010/resequencing/soft/lib';
use lib 'E:/BGI/toGit/perlib/etc';
use strict;
use warnings;
use DBI;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;
#use Galaxy::Data;
#use Fcntl qw(:DEFAULT :flock);

$main::VERSION=0.1.1;

our $opts='i:o:s:d:r:bvf';
our($opt_i, $opt_o, $opt_s, $opt_r, $opt_d, $opt_v, $opt_f, $opt_b);

our $help=<<EOH;
\t-i Input Indel files path (./indel) [ *.gz/*.bz2 files are OK ]
\t-f Input SV files in single file.
\t-s Specie name for shareing SQLite data file (human)
\t  [Specie name MUST be ths SAME throughout !]
\t-d GFF SQLite data file (_dbGFF.sqlite)
\t-o Output SQLite data file (_result_indel.sqlite)
\t-r Output Result txt file (_result_indel.txt)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH
# \t-d dbSNP SQLite data file (_tdbSNP.sqlite)

ShowHelp();

$opt_i='./indel' if ! defined $opt_i;
$opt_d='_dbGFF.sqlite' if ! $opt_d;
$opt_r='_result_indel.txt' if ! $opt_r;
$opt_o='_result_indel.sqlite' if ! $opt_o;
$opt_s='human' if ! $opt_s;

print STDERR "From [$opt_i] [$opt_d]\n  to [$opt_o] [$opt_r]\tSpecie:[$opt_s]\n";
if ($opt_f) {print STDERR "File Input Mode.\n"}
 else {print STDERR "Dir. Input Mode.\n"}
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my @files;
if ($opt_f) {@files=($opt_i);}
 else {
	opendir INDIR,$opt_i  or die "Error: $!\n";
	@files = grep /\.gz$|\.indel\.final$|\.bzip2$/i,readdir INDIR;
	closedir INDIR;
	}

my $start_time = [gettimeofday];
#BEGIN
my %attr = (
    RaiseError => 0,
    PrintError => 1,
    AutoCommit => 0
);

my $dbh = DBI->connect('dbi:SQLite:dbname='.$opt_o,'','',\%attr)
 or die $DBI::errstr;
my $dbhd = DBI->connect('dbi:SQLite:dbname='.$opt_d,'','',\%attr)
 or die $DBI::errstr;

our $sthd=$dbhd->prepare( "SELECT DISTINCT primary_inf,name,start,end FROM gff$opt_s WHERE
 (chrid = :1 AND start <= :2 AND end >= :2) OR
 (chrid = :1 AND start <= :3 AND end >= :3) OR
 (chrid = :1 AND start <= :2 AND end >= :3) OR
 (chrid = :1 AND start >= :2 AND end <= :3)" );

my (%Pinf,%Ainf);

my $sql=q/
CREATE TABLE IF NOT EXISTS reindel{---}
(  chrid TEXT,
   pos INTEGER,
   delta TEXT,
   het TEXT,
   primary_inf TEXT,
   anno_name TEXT );
/;
for (split /;/,$sql) {
	next if /^\s*$/;
	s/{---}/$opt_s/g;
#print "[$_]\n";
	$dbh->do($_) or die $dbh->errstr;
}
$dbh->commit;

my $sth = $dbh->prepare( "INSERT INTO reindel$opt_s
 ( chrid,pos,delta,het,primary_inf,anno_name ) VALUES ( ?,?,?,?,?,? )" );

sub check_mRNA_is_dup($$) {
	my ($mRNA_arr,$elm_arr)=@_;
	my ($mRNA_start,$mRNA_end,$mRNA_name)=@$mRNA_arr[2,3,1];	# primary_inf,name,start,end
	for my $item (@$elm_arr) {
		if ( $$item[1] eq $mRNA_name ) {
			warn "[GFF][!]: element($$item[2]-$$item[3]) outsides mRNA:$mRNA_name($mRNA_start-$mRNA_end)"
			 if ( $mRNA_start > $$item[2] or $mRNA_end < $$item[3] );
			return 1;
		}
	}
	return 0;
}

sub readindel($$) {
	my ($file,$sth)=@_;
	my ($rv,$res,$restr);
	while (<$file>) {
#chromosome_1    241059  I2      AC      *       hete    28      6       10
		my ($chrid,$pos,$delta,$het)=(split /\t/)[0,1,2,5];
		$chrid =~ s/^chr
			(?>
				((?<=^chr)o)?
				((?<=^chro)m)?
				((?<=^chrom)o)?
				((?<=^chromo)s)?
				((?<=^chromos)o)?
				((?<=^chromoso)m)?
				((?<=^chromosom)e)?
			)//xi;
		$delta =~ /^(\D+?)(\d+)$/;
		my ($type,$count,$endpos) = ($1,$2,$pos);
		if ($type =~ /I/i) { $endpos=$pos+$count; }	# well, I/D/etc.
		 elsif ($type =~ /D/i) { $endpos=$pos-$count; }
		my ($start,$end)=sort {$a <=> $b} ($pos,$endpos);
#		($start,$end) = sort {$a <=> $b} ($start,$end);
#warn "$chrid,$svtype,$start,$end\n" if $start > $end;
#print "$chrid,$svtype,$start,$end\n"; sleep 1;
		$sthd->execute($chrid,$start,$end);
		$rv = $sthd->fetchall_arrayref;
		if ($#$rv == -1) {
#			warn "[GFF][x] No infor. for chr:$chrid:$start-$end\n";
			$sth->execute($chrid,$pos,$delta,$het,'UNKNOWN','N/A');
			next;
		}
		my (@Elements,@mRNA)=();
		for (@$rv) {
			#my ($primary_inf,$name)=@$_;
			if ($$_[0] =~ /mRNA/) {push @mRNA,$_;}
			 else {push @Elements,$_;}
		}
		for my $mRNA_arr (@mRNA) {
			next if ! defined $mRNA_arr;
			if ( check_mRNA_is_dup($mRNA_arr,\@Elements) ) { $mRNA_arr = undef; next;}
		}
		# now, no duped mRNA.
		for $res (@Elements,@mRNA) {
			next if ! defined $res;
			$restr="chr:$chrid\t$pos\t$delta\t$het\t$$res[0]\t$$res[1]\t$$res[2]\t$$res[3]\n";
			print OUTFILE $restr;
			++$Pinf{$chrid}{$$res[0]};
			++$Ainf{$$res[0]};
			$sth->execute($chrid,$pos,$delta,$het,@$res[0,1]);
		}
	}
	close $file;
}

my ($t,%pid,$infile,$outfile)=(0);
open(OUTFILE,">",$opt_r) or die "Error: $!\n";

for (@files) {
#	next unless /chr\w{1,2}\./;	# only 22+2+1
	++$t;
	#next if /_random\./;
	$_ = $opt_i.'/'.$_ unless $opt_f;
	print STDERR "parsing [$_] ...\t{$t}\n";
	if (/.bz2$/) {
		open( $infile,"-|","bzip2 -dc $_") or die "Error: $!\n";
	} elsif (/.gz$/) {
    		open( $infile,"-|","gzip -dc $_") or die "Error: $!\n";
	} else {open( $infile,"<",$_) or die "Error: $!\n";}
	readindel($infile,$sth);	# close file within.
}
$dbh->commit;
=pod
print OUTFILE "\n__END__\nSummary:\n";
for my $chr (sort keys %Pinf) {	# {$a <=> $b}
	print OUTFILE "chr:$chr\t";
	for (sort keys %{$Pinf{$chr}}) {
		print OUTFILE $_,":",$Pinf{$chr}{$_},"\t";
	}
	print OUTFILE "\n";
}
=cut
print OUTFILE "\n__END__\nSummary:\n#Chr\t";
my $out;
print OUTFILE "$_\t" for (sort keys %Ainf);
print OUTFILE "\n";
for my $chr (sort keys %Pinf) {	# {$a <=> $b}
	print OUTFILE "$chr\t";
	for (sort keys %Ainf) {
		if (defined $Pinf{$chr}{$_}) {$out=$Pinf{$chr}{$_};}
		 else {$out=0;}
		print OUTFILE $out,"\t";
	}
	print OUTFILE "\n";
}
close OUTFILE;

$sql=q/
CREATE INDEX IF NOT EXISTS svcse{---} ON reindel{---}(chrid,pos);
CREATE INDEX IF NOT EXISTS svpi{---} ON reindel{---}(primary_inf);
CREATE INDEX IF NOT EXISTS svt{---} ON reindel{---}(het);
/;	# chrid,pos,delta,het,primary_inf,anno_name
for (split /;/,$sql) {
	next if /^\s*$/;
	s/{---}/$opt_s/g;
#print "[$_]\n";
	$dbh->do($_) or warn $dbh->errstr;
}
$dbh->commit;

$dbh->disconnect;

$dbhd->rollback;
$dbhd->disconnect;

#my @threads = threads->list();
#$_->join() for @threads;

my $stop_time = [gettimeofday];
#close $infile;
#$|=1;
print STDERR "\n Time Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s)\n";

print STDERR "\033[32;1m Please use [$opt_s] as Specie name in later steps.\033[0;0m\n";
__END__
