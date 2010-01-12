#!/usr/bin/perl -w
#use threads 1.73;
use strict;
use warnings;
use DBI;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;
#use Galaxy::Data;
#use Fcntl qw(:DEFAULT :flock);

$main::VERSION=0.1.1;

our $opts='i:o:s:d:r:bvfk:p:m:t:n';
our($opt_i, $opt_o, $opt_s, $opt_r, $opt_d, $opt_v, $opt_f, $opt_b, $opt_k, $opt_p, $opt_m, $opt_t, $opt_n);

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
\t-k KOG (./Glyma1.KOG)
\t-p PANTHER (./Glyma.PANTHER)
\t-m PFAM (./Glyma.PFAM)
\t-t Samples list (./samples.list) [sample\\s+]
\t-n moNo mode
EOH
# \t-d dbSNP SQLite data file (_tdbSNP.sqlite)

ShowHelp();

$opt_i='./indel' if ! defined $opt_i;
$opt_d='_dbGFF.sqlite' if ! $opt_d;
$opt_r='_result_indel.txt' if ! $opt_r;
$opt_o='_result_indel.sqlite' if ! $opt_o;
$opt_s='human' if ! $opt_s;
$opt_k='./Glyma1.KOG' if ! $opt_k;
$opt_p='./Glyma.PANTHER' if ! $opt_p;
$opt_m='./Glyma.PFAM' if ! $opt_m;
$opt_t='./samples.list' if ! $opt_t;

print STDERR "From [$opt_i] [$opt_d], $opt_k,$opt_p,$opt_m,$opt_t\n  to [$opt_o] [$opt_r]\tSpecie:[$opt_s]\n";
if ($opt_n) {print STDERR "Mono Mode.\n"}
if ($opt_f) {print STDERR "File Input Mode.\n"}
 else {print STDERR "Dir. Input Mode.\n"}
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my (@Samples,$Samples);
open( SAMP,'<',$opt_t) or die "Error: $!\n";
while (<SAMP>) {
	chomp;
	my ($sample)=split /\s+/;
	push @Samples,$sample;
}
close SAMP;
$Samples=join '^',@Samples;

my %hK=();
my %hPa=();
my %hPf=();
getHash($opt_k,\%hK);
getHash($opt_p,\%hPa);
getHash($opt_m,\%hPf);
sub getHash {
	my( $file ,$hash) = @_;
	open (F,"$file") or die $!;
	while (<F>) {
		chomp;
		my($key,$a,$b)=split /\s+/,$_,3;
		$key=~s/\.\d+$//;
		my $value = "#$a\t$b";
	#print "$key\t$value\n";
		$$hash{$key} = $value;
	}
	close F;
}

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
 chrid = :1 AND start <= :2 AND end >= :2" );

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
		my $name=$$item[1];
		if ($$item[0] =~ /gene/i) {
			$mRNA_name=~s/\.\d+$//;
			$name=~s/\.\d+$//;
		}
		if ( $name eq $mRNA_name ) {
			warn "[GFF][!]: element($$item[2]-$$item[3]) outsides mRNA:$mRNA_name($mRNA_start-$mRNA_end)"
			 if ( ($mRNA_start > $$item[2] or $mRNA_end < $$item[3]) and $$item[0] !~ /gene/i);
			return 1;
		}
	}
	return 0;
}

sub readindel($$) {
	my ($file,$sth)=@_;
	my ($rv,$res,$restr);
	while (<$file>) {
		chomp;
		my $line=$_;
		$line =~ s/\^$//;
		my ($chrid,$pos)=split /\t/,$line;
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
		#$delta =~ /^(\D+?)(\d+)$/;
		#my ($type,$count,$endpos) = ($1,$2,$pos);
		#if ($type =~ /I/i) { $endpos=$pos+$count; }	# well, I/D/etc.
		# elsif ($type =~ /D/i) { $endpos=$pos-$count; }
		#my ($start,$end)=sort {$a <=> $b} ($pos,$endpos);
#		($start,$end) = sort {$a <=> $b} ($start,$end);
#warn "$chrid,$svtype,$start,$end\n" if $start > $end;
#print "$chrid,$svtype,$start,$end\n"; sleep 1;
		$sthd->execute($chrid,$pos);
		$rv = $sthd->fetchall_arrayref;
		if ($#$rv == -1) {
#			warn "[GFF][x] No infor. for chr:$chrid:$start-$end\n";
			$sth->execute($chrid,$pos,'','','UNKNOWN','N/A');
			my $tt=".\t.";
			print OUTFILE "$line\tUnknown\tN/A\t$tt\t$tt\t$tt\n";
			next;
		}
		my (@Elements,@mRNA,@EXON,@CDS,@UTR,@Gene,@Intron)=();	# $opt_n
		for (@$rv) {
			#my ($primary_inf,$name)=@$_;
			if ($$_[0] =~ /mRNA/i) {push @mRNA,$_;}
			if ($$_[0] =~ /exon/i) {push @EXON,$_;}
			if ($$_[0] =~ /CDS/i) {push @CDS,$_;push @Elements,$_;}
			if ($$_[0] =~ /UTR/i) {push @UTR,$_;push @Elements,$_;}
			if ($$_[0] =~ /Gene/i) {push @Gene,$_;}
			if ($$_[0] =~ /Intron/i) {push @Intron,$_;push @Elements,$_;}
			 else {push @Elements,$_;}
		}
		unless ($opt_n) {	# NOT WORKING !!!
		for my $mRNA_arr (@mRNA) {
			next if ! defined $mRNA_arr;
			if ( check_mRNA_is_dup($mRNA_arr,\@Elements) ) { $mRNA_arr = undef; next;}
		}
		for my $mRNA_arr (@Gene) {
			next if ! defined $mRNA_arr;
			if ( check_mRNA_is_dup($mRNA_arr,\@Elements) ) { $mRNA_arr = undef; next;}
		}
		for my $mRNA_arr (@EXON) {
			next if ! defined $mRNA_arr;
			if ( check_mRNA_is_dup($mRNA_arr,\@Elements) ) { $mRNA_arr = undef; next;}
		}
		}
		# now, no duped mRNA.
		for $res (@CDS,@Intron,@UTR,@EXON,@mRNA,@Gene) {
			next if ! defined $res;
			my $key=$$res[1];
			$key =~s/\.\d+$//;
			my ($kk,$pp,$mm);
			$kk=(defined $hK{$key})?$hK{$key}:".\t.";
			$pp=(defined $hPa{$key})?$hPa{$key}:".\t.";
			$mm=(defined $hPf{$key})?$hPf{$key}:".\t.";
			$restr="$line\t$$res[0]\t$$res[1]\t$kk\t$pp\t$mm\n";
			print OUTFILE $restr;
			++$Pinf{$chrid}{$$res[0]};
			++$Ainf{$$res[0]};
			$sth->execute($chrid,$pos,'','',@$res[0,1]);
			last if $opt_n;
		}
	}
	close $file;
}

my ($t,%pid,$infile,$outfile)=(0);
open(OUTFILE,">",$opt_r) or die "Error: $!\n";
print OUTFILE "ChrID\tPos\tUp200\tDown200\t$Samples\tType\tGene\tKOG_ID\tKOG_nfo\tPANTHER_ID\tPANTHER_nfo\tPFAM_ID\tPFAM_nfo\n";

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
=pod
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
=cut
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
