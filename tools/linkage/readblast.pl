#!/bin/env perl
use strict;
use warnings;

unless (@ARGV > 0) {
    print "perl $0 <max trim> <outfile> <infiles>\n";
    exit 0;
}
# infiles should be of :
# -evalue 3e-11 -best_hit_overhang 0.25 -best_hit_score_edge 0.05 -outfmt '6 qseqid sseqid pident length nident mismatch gapopen qstart qend sstart send evalue bitscore btop'
# since each marker can only in one blast result file, we can split by chr.
# For PA64 g=347534994(non-N), 0.01*1/g=2.877e-11. Alpha=0.96 %

my $maxTrim=shift;
my $out_file = shift;
my @Sinf;

sub getO($$$) {
	my ($col,$m,$arr)=@_;
	#my @A=sort {$a->[$col] <=> $b->[$col]} @$arr;
	my %A;
	push @{$A{ $_->[$col] }},$_ for @$arr;
	my $id=(sort {$a <=> $b} keys %A)[$m];	# ASC
	return $A{$id};
}
sub pasteSinf() {
	#my ()=@_;
	for (@Sinf) {
		my ($Qid,$Sid,$Pidentity,$AlnLen,$identical,$mismatches,$Gap,$Qs,$Qe,$Ss,$Se,$E,$bitScore,$BTOP)=@$_;
		my @t=split //,$BTOP;
		for (0..$#t-1) {
			++$identical if $t[$_] =~ /[MRWSYK]/i and $t[$_+1] ne '-';
		}
		$_=[$Qid,$Sid,$Pidentity,$AlnLen,$identical,$mismatches,$Gap,$Qs,$Qe,$Ss,$Se,$E,$bitScore,$BTOP];
	}
	my $A=&getO(4,-1,\@Sinf);	# by col.4 identical, max
	$A=&getO(11,0,$A);	# by col.11 E, min
	return [] if $#$A > 0;	# if still more than 1, is duplicated
	return $$A[0];
}

if (-s $out_file) {
	die "[x]$out_file exists, remove first !\n";
} else {
	warn "[!]Open $out_file for write.\n";
}
open O,'>',$out_file or die $!;
open LOG,'>',$out_file.'.log' or die $!;
while (my $in_file=shift) {
	print LOG ">$in_file\t";
	open I,'<',$in_file or warn "\n[!]Error opening $in_file: $!\n";
	my ($lastQid,$MarkerLen);
	while (<I>) {
		next if /^#/;
# Fields: query id, subject id, % identity, alignment length, identical, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, BTOP
#-outfmt '6 qseqid sseqid pident length nident mismatch gapopen qstart qend sstart send evalue bitscore btop'
		chomp;
		my ($Qid,$Sid,$Pidentity,$AlnLen,$identical,$mismatches,$Gap,$Qs,$Qe,$Ss,$Se,$E,$bitScore,$BTOP)=split /\t/;
		my @Dat=($Qid,$Sid,$Pidentity,$AlnLen,$identical,$mismatches,$Gap,$Qs,$Qe,$Ss,$Se,$E,$bitScore,$BTOP);
		if ($lastQid and ($lastQid ne $Qid)) {
			#print STDERR "-";
			my @t=@{&pasteSinf()};
			goto tNEXT if $#t == -1;
			my ($Qs,$Qe)=@t[7,8];	# @{$t[0]}[7,8]
			if ( ($Qs > $maxTrim) or (($MarkerLen-$Qe)> $maxTrim) ) {
				print LOG "|$Qs,$Qe\n";
				goto tNEXT;
			} else {
				print LOG "-";
			}
			print O join("\t",@t,scalar @Sinf),"\n";
tNEXT:		@Sinf=(\@Dat);
		} else {
			push @Sinf,\@Dat;
		}
		$lastQid=$Qid;
		$MarkerLen=1 + 2*(split /m/,$Qid)[-1];
	}
	close I;
	print LOG ".\n";
}
close LOG;
close O;

__END__
./readblast.pl 15 t.out ex32Chr01.pa64m7e
cat chrorder|while read a;do ./readblast.pl 15 ./out/f45$a.pa64 ./out/ex45$a.pa64m6;done	# 2> ./out/f45$a.pa64.log

perl -lane 'BEGIN {my %d} ++$d{$F[1]}; END {print "$_\t$d{$_}" for sort {$d{$a} <=> $d{$b}} keys %d;}' ./out/f45*.pa64 | tee stat.pa64
grep Scaffold stat.pa64 |perl -lane 'BEGIN {my %a;my $b} ++$a{$F[1]};$b += $F[1]; END {print "$_\t$a{$_}" for sort {$a<=>$b} keys %a; print "\n$b"}'

#!/bin/sh
#$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH,BLASTDB
#$ -cwd -r y -l vf=900m
#$ -o /dev/null -e /dev/null
#$ -S /bin/bash -t 1-12
SEEDFILE=./doblast.cmd
SEED=$(sed -n -e "$SGE_TASK_ID p" $SEEDFILE)
export BLASTDB=/ifs1/POPULATION/Rice/gapfilling/blast
eval ./blastn -task megablast -db PA64scaf -outfmt \'6 qseqid sseqid pident length nident mismatch gapopen qstart qend sstart send evalue bitscore btop\' -num_threads 2 -evalue 3e-11 -best_hit_overhang 0.25 -best_hit_score_edge 0.05 $SEED

$ cat doblast.cmd
-query ./marker/ex45Chr01.marker -out ./out/ex45Chr01.pa64m6
-query ./marker/ex45Chr02.marker -out ./out/ex45Chr02.pa64m6
