#!/bin/env perl
use lib '/share/raid010/resequencing/soft/lib';
use lib 'E:/BGI/toGit/perlib/etc';
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;
use Bio::SeqIO;
use Bio::Seq::Quality;

$main::VERSION=0.0.1;


our $opts='i:q:o:s:l:r:b';
our($opt_i, $opt_q, $opt_s, $opt_l, $opt_r, $opt_o, $opt_b);

#our $desc='1.filter fq, 2.stats';
our $help=<<EOH;
\t-i Input FASTA file
\t-q Quality file
\t-o Output Prefix for (./out)_{1,2}.fq
\t-s Insert_Size (500)
\t-r Reads_Len (75)
\t-l slip_dis (50)
\t-b No pause for batch runs
EOH

ShowHelp();

die "[x]Must specify FASTA file !\n" if ! $opt_i;
die "[x]-i $opt_i not exists !\n" unless -f $opt_i;
$opt_s = 500 unless $opt_s;
$opt_r = 75 unless $opt_r;
$opt_l = 50 unless $opt_l;
$opt_o = './out' unless $opt_o;

my $min_r=33;

print STDERR "From [$opt_i],[$opt_q] with [$opt_s,$opt_r,$opt_l] to [$opt_o]_{1,2}.fq\n";
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
my ($fafile,$qfile);
if ($opt_i =~ /\.gz$/) {
	$fafile = "gzip -dc $opt_i |";
} elsif ($opt_i =~ /\.bz2$/) {
 	$fafile = "bzip2 -dc $opt_i |";
} else { $fafile=$opt_i; }

if ($opt_q =~ /\.gz$/) {
	$qfile = "gzip -dc $opt_q |";
} elsif ($opt_q =~ /\.bz2$/) {
 	$qfile = "bzip2 -dc $opt_q |";
} else { $qfile=$opt_q; }

my $in_seq_obj =
  Bio::SeqIO->new( -file => $fafile,
		   -format => 'fasta',
		 );

my $in_qual_obj =
  Bio::SeqIO->new( -file => $qfile,
		   -format => 'qual',
		 );

open FQ1,'>',$opt_o.'_1.fq' or die "[x]Error on ${opt_o}_{1,2}.fq: $!\n";
open FQ2,'>',$opt_o.'_2.fq';

sub revcom($) {
	my $str = $_[0];
	$str =~ tr/acgtrymkswhbvdnxACGTRYMKSWHBVDNX/tgcayrkmswdvbhnxTGCAYRKMSWDVBHNX/;
	my $rev = reverse $str;
	#$rev =~ tr/[](){}<>/][)(}{></;
	return $rev;
}
sub Qstr($) {
	my @arr=@{$_[0]};
	my (@res,$Q);
#$Q = 10 * log( 1 + 10 ** (ord($sq) - 64) / 10.0) ) / log(10);
	for my $q (@arr) {
		#$e = 10 ** (-$q/10);
		#$sQ = -10 * log($e / (1 - $e)) / log(10);
		$Q=chr($q+64);	# version 1.3+ of the GA pipeline encodes Phred quality scores from 0-62 using ASCII 64-126
		push @res,$Q;
	}
	return join('',@res);
}

print STDERR '[!]Begin: ';
my ($Count,$as,$bs,$aq,$bq,$t)=(0);
while (1) {
	my $seq_obj  = $in_seq_obj->next_seq || last;
	my $qual_obj = $in_qual_obj->next_seq;
	die "[x]foo!\n" unless $seq_obj->id eq $qual_obj->id;
	my $id = $seq_obj->id;
	my $seq = $seq_obj->seq;
	my $qual = &Qstr($qual_obj->qual);
	my $len=length $seq;
	next if $len < $min_r * 2;
	if ($len <= 2*$opt_r) {
		$t=int(0.5+$len/2);
		$as=substr $seq,0,$t;
		$bs=&revcom(substr $seq,$t-1);
		$aq=substr $qual,0,$t;
		$bq=reverse(substr $qual,$t-1);
		print FQ1 "\@${id}_${len}/1\n$as\n+\n$aq\n";
		print FQ2 "\@${id}_${len}/2\n$bs\n+\n$bq\n";
	} else {	# $len mod $opt_l is another border, skipped
		my $s=0;
		while ($len-$s > $opt_s) {
			$as=substr $seq,$s,$opt_r;
			$bs=&revcom(substr $seq , $s+$opt_s-$opt_r , $opt_r);
			$aq=substr $qual,$s,$opt_r;
			$bq=reverse(substr $qual , $s+$opt_s-$opt_r , $opt_r);
			if (length $bs < 33) {
				$bs .= 'N' x 32;
				$bq .= '~' x 32;
			}
			print FQ1 "\@${id}_${len}_${s}/1\n$as\n+\n$aq\n";
			print FQ2 "\@${id}_${len}_${s}/2\n$bs\n+\n$bq\n";
			$s += $opt_l;
		}
	}
	++$Count;
	print STDERR '.' unless $Count % 1000;
}
warn "\n[!] Total $Count PE out.\n";
close FQ2;
close FQ1;


#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__
find ./oryza_sativa__indica_cultivar_group/fasta.oryza_sativa__indica_cultivar_group.*.gz|perl -lane '(undef,$p,$f)=split /\//;@t=split /\./,$f;print "-i ./$p/$f -q ./$p/qual.$t[1].$t[2].gz -o ./fq/fa9311_$t[2] 2>./fq/fa9311_$t[2].log"' > cmd.lst

#!/bin/sh
#$ -N "EmuFQ"
#$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#$ -cwd -r y -l vf=116M
#$ -o /dev/null -e /dev/null
#$ -S /bin/bash -t 1-13
SEEDFILE=./cmd.lst
SEED=$(sed -n -e "$SGE_TASK_ID p" $SEEDFILE)
eval perl ./fa2fq.pl -b $SEED

./soappe.pl 495,505 200 ./fq/fa9311_001_1.fq ./fq/fa9311_001_2.fq ./9311/9311_main_chr.fa.index ./soap/t1

find ./fq/*_1.fq|perl -lane '(undef,$p,$f)=split /\//;@t=split /_1\./,$f;print "./$p/$f ./$p/$t[0]_2.fq ./9311/9311_main_chr.fa.index ./soap/$t[0]"' > cmds.lst

#!/bin/sh
#$ -N "EmuSoap"
#$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#$ -cwd -r y -l vf=1.8g
#$ -o /dev/null -e /dev/null
#$ -hold_jid EmuFQ
#$ -S /bin/bash -t 1-13
SEEDFILE=./cmds.lst
SEED=$(sed -n -e "$SGE_TASK_ID p" $SEEDFILE)
eval perl ./soappe.pl 495,505 200 $SEED


use Bio::SeqIO;
use Bio::Seq::Quality;

use Getopt::Long;

die "pass a fasta and a fasta-quality file\n"
  unless @ARGV;


my ($seq_infile,$qual_infile)
  = (scalar @ARGV == 1) ?($ARGV[0], "$ARGV[0].qual") : @ARGV;

## Create input objects for both a seq (fasta) and qual file

my $in_seq_obj =
  Bio::SeqIO->new( -file   => $seq_infile,
		   -format => 'fasta',
		 );

my $in_qual_obj =
  Bio::SeqIO->new( -file   => $qual_infile,
		   -format => 'qual',
		 );

my $out_fastq_obj =
  Bio::SeqIO->new( -format => 'fastq'
		 );

while (1){
  ## create objects for both a seq and its associated qual
  my $seq_obj  = $in_seq_obj->next_seq || last;
  my $qual_obj = $in_qual_obj->next_seq;

  die "foo!\n"
    unless
      $seq_obj->id eq
	$qual_obj->id;

  ## Here we use seq and qual object methods feed info for new BSQ
  ## object.
  my $bsq_obj =
    Bio::Seq::Quality->
	new( -id   => $seq_obj->id,
	     -seq  => $seq_obj->seq,
	     -qual => $qual_obj->qual,
	   );

  ## and print it out.
  $out_fastq_obj->write_fastq($bsq_obj);
}

     my $seqin = Bio::SeqIO->new(-file   => "/usr/local/bin/gunzip -c $infile |",
                                 -format => $informat);

     my $seqout = Bio::SeqIO->new(-file   => ">$outfile",
                                  -format => 'Fasta');

