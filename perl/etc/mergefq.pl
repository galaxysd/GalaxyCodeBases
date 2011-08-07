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
	#chomp(my $a=<$fh>) or last;
	chomp(my $genome=<$fh>) or last;
	chomp(my $c=<$fh>) or last;
	chomp(my $d=<$fh>) or last;
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
find . -name '1.fq.gz.out'|perl -lane '$n=$_;s/1\.fq\.gz\.out$//;$p=$_;print "./bwa aln -n $_ -o 0 -N -I ${p}1merge.fa $n >${p}1.$_.sai 2>${p}1.$_.sai.log" for (3,6,10,2,5,1,20)' > rsai.cmd
find . -name '1.fq.gz.out'|perl -lane '$n=$_;s/1\.fq\.gz\.out$//;$p=$_;print "./bwa samse -n 500 -f ${p}1.$_.sam ${p}1merge.fa ${p}1.$_.sai $n 2>${p}1.$_.sam.log" for (3,6,10,2,5,1,20)' > rsam.cmd

$ cat rdo1.sh
#!/bin/sh
#$ -N "snp_1"
#$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#$ -cwd -r y -l vf=500m,p=1
#$ -o /dev/null -e /dev/null
#$ -S /bin/bash -t 1-21
SEEDFILE=./rsai.cmd
SEED=$(sed -n -e "$SGE_TASK_ID p" $SEEDFILE)
WCNS=./rsam.cmd
STR=$(sed -n -e "$SGE_TASK_ID p" $WCNS)
SAMFILE=$(echo $STR|awk '{print $6}')

eval $SEED
eval $STR
grep -P '\tXT:A:R\t' $SAMFILE >$SAMFILE.r
