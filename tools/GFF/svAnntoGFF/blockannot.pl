#!/usr/bin/perl -w
#use threads 1.73;
use strict;
use warnings;
use lib '/home/huxuesong/perl/lib64/perl5/5.8.5/x86_64-linux-thread-multi/';
use lib '/home/huxuesong/perl/lib64/perl5/site_perl/5.8.5/x86_64-linux-thread-multi/';
use lib '/home/huxuesong/gentoo/home/PERL5LIB/';
use DBI;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;
#use Galaxy::Data;
#use Fcntl qw(:DEFAULT :flock);

$main::VERSION=0.1.1;

our $opts='i:o:s:d:r:bvf';
our($opt_i, $opt_o, $opt_s, $opt_d, $opt_v, $opt_f, $opt_b);

our $help=<<EOH;
\t-i Input block files path (./ex)
\t-f Input block files in single file.
\t-s Specie name for shareing SQLite data file (sorghum)
\t  [Specie name MUST be ths SAME throughout !]
\t-d GFF SQLite data file (_dbGFF.sqlite)
\t-o Output file suffix (.annot)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH
# \t-d dbSNP SQLite data file (_tdbSNP.sqlite)

ShowHelp();

$opt_i='./ex' if ! defined $opt_i;
$opt_d='_dbGFF.sqlite' if ! $opt_d;
$opt_o='.annot' if ! $opt_o;
$opt_s='sorghum' if ! $opt_s;

print STDERR "From [$opt_i] [$opt_d]\n  to [$opt_o]\tSpecie:[$opt_s]\n";
if ($opt_f) {print STDERR "File Input Mode.\n"}
 else {print STDERR "Dir. Input Mode.\n"}
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my @files;
if ($opt_f) {@files=($opt_i);}
 else {
	opendir INDIR,$opt_i  or die "Error: $!\n";
	@files = grep /\.block$/i,readdir INDIR;
	closedir INDIR;
	}

my $start_time = [gettimeofday];
#BEGIN
my %attr = (
    RaiseError => 0,
    PrintError => 1,
    AutoCommit => 0
);

my $dbhd = DBI->connect('dbi:SQLite:dbname='.$opt_d,'','',\%attr)
 or die $DBI::errstr;

our $sthd=$dbhd->prepare( "SELECT DISTINCT primary_inf,name FROM gff$opt_s WHERE
 (chrid = :1 AND start <= :2 AND end >= :2)" );

#my (%Pinf,%Ainf);

sub check_mRNA_is_dup($$) {
	my ($mRNA_arr,$elm_arr)=@_;
	my ($mRNA_name)=@$mRNA_arr[1];	# primary_inf,name,start,end
	for my $item (@$elm_arr) {
		if ( $$item[1] eq $mRNA_name ) {
			return 1;
		}
	}
	return 0;
}

sub readSV($) {
	my ($file)=@_;
	my ($rv,$res,$restr);
	while (<$file>) {
#chromosome_5    1394    1395    2
		my ($chrid,$start,$end,$len)=split /\t/;
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
		($start,$end) = sort {$a <=> $b} ($start,$end);	# EP ?
#warn "$chrid,$svtype,$start,$end\n" if $start > $end;
#print "$chrid,$svtype,$start,$end\n"; sleep 1;
		for my $pos ($start..$end) {
			$sthd->execute($chrid,$pos);
			$rv = $sthd->fetchall_arrayref;
			if ($#$rv == -1) {
#				warn "[GFF][x] No infor. for chr:$chrid:$start-$end\n";
				print OUTFILE "chr:$chrid\t$pos\tnongene\tNA\n";
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
				$restr="chr:$chrid\t$pos\t$$res[0]\t$$res[1]\n";
				print OUTFILE $restr;
#				++$Pinf{$chrid}{$$res[0]};
#				++$Ainf{$$res[0]};
			}
		}
	}
	close $file;
}
my ($t,%pid,$infile,$outfile)=(0);

for (@files) {
#	next unless /chr\w{1,2}\./;	# only 22+2+1
	++$t;
	#next if /_random\./;
	$_ = $opt_i.'/'.$_ unless $opt_f;
	print STDERR "parsing [$_] ...\t{$t}\n";
	open( $infile,"<",$_) or die "Error: $!\n";
	open(OUTFILE,">",$_.$opt_o) or die "Error: $!\n";
	readSV($infile);	# close file within.
	close OUTFILE;
}

$dbhd->rollback;
$dbhd->disconnect;

my $stop_time = [gettimeofday];
#close $infile;
#$|=1;
print STDERR "\n Time Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s)\n";

print STDERR "\033[32;1m Please use [$opt_s] as Specie name in later steps.\033[0;0m\n";
__END__
