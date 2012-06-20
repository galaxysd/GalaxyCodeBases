#!/bin/env perl
use strict;
use warnings;

die "Usage: $0 <qusb_job.sh>\n" if @ARGV != 1;
my ($qsh)=@ARGV;

die "sh name [$qsh] cannot begin with '.' !\n" if $qsh =~ /^\./;
die "sh name [$qsh] must be plain filename !\n" if $qsh =~ /[\/\\]/;
die "sh name [$qsh] contains no words !\n" unless $qsh =~ /\w/;

print STDERR "From [$qsh] to [.$qsh/*.sh], output path will be DELETED first.\n";
print STDERR 'press [Enter] to continue...'; <STDIN>;

open IN,'<',$qsh or die "Cannot open [$qsh]: $!\n";
system('rm','-vfr',".$qsh");
system('mkdir',".$qsh");
my (@Heads,@Lines);
my $flag = 1;
my ($jaBegin,$jaEnd,$jaLine)=(0,0,'');
while (<IN>) {
	chomp;
	if ($flag && /^#[!\$]/) {
		if (my $ret = search_jobarray($_)) {
			($jaBegin,$jaEnd)=@$ret;
			$jaLine = $_;
			s/(\b| )-t (\d+)-(\d+)\b//;
		}
		push @Heads,$_;
		die "Only 1 \"-t\" supported !\n" if /(\b| )-t /;
	} else {
		$flag = 0;
		die "Cannot set SGE_TASK_ID in [$_]\n" if (! /^\s*#/) and /SGE_TASK_ID=/;
		push @Lines,$_;
	}
}
close IN;

print '=' x 20,' HEADER ','=' x 21,"\n";
print '[',join("]\n[",@Heads),"]\n";
print '=' x 20,' CONTENT ','=' x 20,"\n";
print '[',join("]\n[",@Lines),"]\n";
print '=' x 49,"\n";
if ($jaBegin) {
	print "ArrayJob: [$jaBegin] - [$jaEnd] found in [$jaLine]!\n";
} else {
	print "Not an ArrayJob !\nqsub directly.\n";
	exit;
}

for my $ajobid ($jaBegin .. $jaEnd) {
	open OUT,'>',".$qsh/j${ajobid}_$qsh" or die "$!"; # a valid object name (cannot start with a digit).
	print OUT join("\n",@Heads);
	print OUT "\n\nSGE_TASK_ID=$ajobid\n\n";
	print OUT join("\n",@Lines),"\n";
	close OUT;
}
print "Run: find \".$qsh/\" -name \"*_$qsh\"|while read a;do qsub \"\$a\";done\n ";

sub search_jobarray() {
	my $line = $_[0];
	if ($line =~ /(\b| )-t (\d+)-(\d+)/) {
		return [$2,$3];
	} else {
		return 0; # NULL
	}
}

__END__

#!/bin/bash
#$ -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH
#$ -cwd -r y -l vf=13g,p=1
#$ -o /dev/null -e /dev/null
#$ -S /bin/bash -t 1-2
# c test
 # SGE_TASK_ID=1
SEEDFILE=sam.lst
SEED=$(sed -n -e "${SGE_TASK_ID} p" $SEEDFILE)
MAIN=`basename "$SEED"`
# comments here
./pIRS/baseCalling_Matrix_calculator -l 101 -r ./ref/hg19.fa -s dbsnp132with1kgenome.lst.bz2 -o S_${SGE_TASK_ID} -c chr22 -b $SEED > s_${SGE_TASK_ID}.log 2> s_${SGE_TASK_ID}.err
#ccc

