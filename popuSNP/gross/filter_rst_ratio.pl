#!/usr/bin/perl -w
use strict;

print STDERR $0,"\n";

my $indir = shift;
my $chr = shift;
$indir =~ s/\/$//g;

my $rst = "$indir/$chr\.rst";
my $ratio = "$indir/$chr\.ratioCheck";

my %hashrst;
open (A,$rst) || die $!;
while(<A>){
	chomp;
	my @line=split;
	if ($line[2] >= 0.005 && $line[4] >= 0){
		$hashrst{"$line[0]\t$line[1]"}=1;
	}
}
close A;

open (A,$ratio) || die $!;
while(<A>){
	chomp;
	my @t = split(/\t/);
	if (exists $hashrst{"$t[0]\t$t[1]"}){
#	if ($t[-1] <= 2 && exists $hashrst{"$t[0]\t$t[1]"}){
		print join ("\t",@t[0..$#t]),"\n";
	}
}
close A;

__END__
cat chrorder | perl -lane 'open O,">./shell/${_}_FR.sh";$a="./population/".$_.".population.snp" ;print O "\#\$ -N FR_$_ -hold_jid RA_$_,LC_$_ -cwd -r y -l vf=1g,p=1 -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH -o /dev/null -e $a.err";print O "./filter_rst_ratio.pl ./population/ $_ > $a";close O;'
