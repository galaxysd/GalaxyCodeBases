#!/bin/env perl
use strict;
use warnings;
#use Data::Dump;

unless (@ARGV){
	print "perl $0 <dat.gz>\n";
	exit;
}

my ($in)=@ARGV;
my $out=$in;
$out =~ s/\.dat\.gz$//;
$out =~ /^(.+\/)([^\/]+)$/;
my ($path,$main)=($1,$2);
warn "From [$main] @ [$path]\n";
open IN,'-|',"gzip -dc $in" or die "[x]Error opening [$in] with gzip: $!\n";
$_=<IN>;
chomp;
s/^#//;
my ($FQa,$FQb)=split /,/;
open FQA,'-|',"gzip -dc $FQa" or die "[x]Error opening [$FQa] with gzip: $!\n";
open FQB,'-|',"gzip -dc $FQb" or die "[x]Error opening [$FQb] with gzip: $!\n";

sub readFQs($) {
	chomp( my ($ida,$seqa,undef,$qa,$idb,$seqb,undef,$qb)=(scalar <FQA>,scalar <FQA>,scalar <FQA>,scalar <FQA>,scalar <FQB>,scalar <FQB>,scalar <FQB>,scalar <FQB>) );
	unless (defined $qa) {
		unless (defined $qb) {
			return -1;
		}
		return -2;
	} elsif (! defined $qb) {
		return -2;
	} else {
		$_[0]=[ [$ida,$seqa,$qa],[$idb,$seqb,$qb] ];
#dd $_[0];
		return 0;
	}
	#$ida=<FQA>; $seqa=<FQA>; <FQA>; $qa=<FQA>;
	#$ida=<FQA>; $seqa=<FQA>; <FQA>; $qa=<FQA>;
}
sub writeFQs($$) {
	my ($Dat,$FHs)=@_;
	my ($a,$b)=@$FHs;
	my ($ida,$seqa,$qa)=@{$$Dat[0]};
	my ($idb,$seqb,$qb)=@{$$Dat[1]};
	print $a "$ida\n$seqa\n+\n$qa\n";
	print $b "$idb\n$seqb\n+\n$qb\n";
#warn "$idb\n$seqb\n+\n$qb\n";
}
sub getOutFHs($) {
	my ($prefix)=@_;
	my ($a,$b);
	open $a,'|-',"gzip -9c - >${path}${prefix}_${main}_1.fq.gz" or die "[x]Error opening [${path}${prefix}_${main}_1.fq.gz] with gzip: $!\n";
	open $b,'|-',"gzip -9c - >${path}${prefix}_${main}_2.fq.gz" or die "[x]Error opening [${path}${prefix}_${main}_2.fq.gz] with gzip: $!\n";
	return [$a,$b];
}
sub closeOutFH($) {
	my ($a,$b)=@{$_[0]};
	close $a;
	close $b;
}

my %FHs=(
	0 => &getOutFHs('U'),
	1 => &getOutFHs('A'),
	2 => &getOutFHs('B'),
	3 => &getOutFHs('H'),
);
#my ($FHsA,$FHsB,$FHsH,$FHsU)=( &getOutFHs('A'),&getOutFHs('B'),&getOutFHs('H'),&getOutFHs('U'), );
<IN>;
my ($lastid,$lasttype)=(0,0);
my $Dat;
while (<IN>) {
	chomp;
	my ($id,$type)=split /\t/;
	if ($lastid < $id) {
		for ($lastid .. $id-1) {
			readFQs($Dat);# or warn "[!]Error: $!";
#dd $Dat;
			writeFQs($Dat,$FHs{$lasttype});
		}
	}
	($lastid,$lasttype)=($id,$type);
}
while (readFQs($Dat)==0) {
	writeFQs($Dat,$FHs{$lasttype});
}
closeOutFH($_) for values %FHs;
close IN;
close FQA;
close FQB;

__END__
#!/bin/sh
#$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#$ -cwd -r y -l vf=600m
#$ -o /dev/null -e /dev/null
#$ -S /bin/bash -t 1-195
SEEDFILE=./dat2.lst
SEED=$(sed -n -e "$SGE_TASK_ID p" $SEEDFILE)
perl ./extractFQ.pl "$SEED"
