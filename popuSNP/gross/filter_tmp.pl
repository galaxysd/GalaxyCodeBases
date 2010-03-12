#!/usr/bin/perl -w
use strict;
print STDERR $0,"\n";

my %mark;

die "perl $0 <population.snp> <chr_pos file(s)>\n" if @ARGV < 2;
print STDERR $0,"\n";

my $infile = shift;
my $t='';

while(<>) {
	if ($ARGV ne $t) {	# for debug only. It is already a silly slow programme.
		warn "[!]< $ARGV\n";
		$ARGV = $t;
	}
	my ($chr,$pos)=split /\t/;
	++$mark{$chr}{$pos};
}
warn "\n[!]>\n";

open (A,$infile) || die $!; ##genotyping result
while(<A>){
	#chomp;
	my ($chr,$pos)=split /\t/;
	next unless defined $mark{$chr}{$pos};
	print;
}
close A;

__END__
cat chrorder | perl -lane 'open O,">./shell/${_}_RA.sh";$a="./population/".$_.".ratioCheck" ;print O "\#\$ -N RA_$_ -hold_jid LC_$_ -cwd -r y -l vf=1g,p=1 -v PERL5LIB,PATH,PYTHONPATH,LD_LIBRARY_PATH -o /dev/null -e $a.err";print O "./check_allSNP.pl ./population/${_}.add_cn ./population/${_}.LC > $a";close O;'
