#!/bin/env perl
use strict;
use warnings;

die "Usage: $0 <mergedfa> <inputfa>\n" if @ARGV != 2;
my ($out,$in)=@ARGV;
warn "From [$in] to [$out]\n";

sub openfile($) {
    my ($filename)=@_;
    my $infile;
    if ($filename=~/.bz2$/) {
	    open( $infile,"-|","bzip2 -dc $filename") or die "Error opening $filename: $!\n";
    } elsif ($filename=~/.gz$/) {
     	open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
    } else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
    return $infile;
}

my $fh=&openfile($in);
my $SEQ='';
my $Ncnt=0;
while (<$fh>) {
    s/^>//;
	/^(\S+)/ or next;
	my $seqname = $1;
    #print STDERR " >$seqname ...";
	$/=">";
	my $genome=<$fh>;
	chomp $genome;
	$genome=~s/\s//g;
	$/="\n";
    $SEQ .= $genome;
    my $n=($genome=~s/[^ATCG]/A/ig);
    $Ncnt += $n;
    #print STDERR "\b\b\b   \t",length $genome,".\n";
	$genome='';
}

open OUT,'>',$out or die "Error opening $out: $!\n";
print OUT ">merged $Ncnt\n$SEQ\n";
close OUT;

__END__
./bwa aln -n 3 -o 0 -I ./ant/ant_meg_free.fa ./ant/ant_error_free_100_500_1.fa >./ant/ant_meg_free.3.sai 2>./ant/ant_meg_free.3.sai.log
./bwa samse -n 5 -f ./ant/ant_meg_free.3.sam ./ant/ant_meg_free.fa ./ant/ant_meg_free.3.sai ./ant/ant_error_free_100_500_1.fa 2>./ant/ant_meg_free.3.sam.log

XT:A:R

XA:Z:merged,+10484501,100M,0;
XA:Z:merged,-15123200,100M,1;

find . -name '*meg*.fa'|perl -lane '/\/(\w+)\//;next if /ant/;system "./bwa index $_"'
find . -name '*_1.fa'|perl -lane '/\/(\w+)\//;$p=$1;$n=$_;s/_error_([^_]+)_.+/_err_$1/;$err=$1;$m=$_;print STDERR "${m}_?";print "./bwa aln -n $_ -o 0 -N -I ./$p/${p}_meg_$err.fa $n >$m.$_.sai 2>$m.$_.sai.log" for (3,6,10,2,5,1,20)' > sai.cmd

find . -name '*_1.fa'|perl -lane '/\/(\w+)\//;$p=$1;$n=$_;s/_error_([^_]+)_.+/_err_$1/;$err=$1;$m=$_;print STDERR "${m}_?";print "./bwa samse -n 500 -f $m.$_.sam ./$p/${p}_meg_$err.fa $m.$_.sai $n 2>$m.$_.sam.log" for (3,6,10,2,5,1,20)' > sam.cmd

SAMFILE=$(echo $SEED|awk '{print $6}')

grep -P '\tXT:A:R\t' $SAMFILE >$SAMFILE.r

./bwa aln -n 6 -o 0 -I ./ant/ant_meg_free.fa ./ant/ant_error_free_100_500_1.fa >./ant/ant_meg_free.6.sai 2>./ant/ant_meg_free.6.sai.log
./bwa aln -n 10 -o 0 -I ./ant/ant_meg_free.fa ./ant/ant_error_free_100_500_1.fa >./ant/ant_meg_free.10.sai 2>./ant/ant_meg_free.10.sai.log

./bwa samse -n 500 -f ./ant/ant_err_1%.3.sam ./ant/ant_meg_1%.fa ./ant/ant_err_1%.3.sai ./ant/ant_error_1%_100_500_1.fa 2>./ant/ant_err_1%.3.sam.log
