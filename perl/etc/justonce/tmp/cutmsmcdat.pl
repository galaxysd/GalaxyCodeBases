#!/usr/bin/perl -w

use strict;

my ($fname) = @ARGV;

$fname =~ /\b(chr\w+)\b/;
my $chrid = $1;

my @indcnt = 1 .. 9;

warn "($1,$fname)\n";
my (@OutFH,@remained);

open IN,'<',$fname or die $!;

for (@indcnt) {
	my $t;
	open $t,'>',"dat/${chrid}.s$_.msmcdat" or die $!;
	$OutFH[$_] = $t;
	$remained[$_] = 0;
}

while (<IN>) {
	chomp;
	my @dat = split /\t/,$_;
	#print "[$dat[-1]]\n";
	next if length $dat[-1] != 54;
	next unless $dat[2] > 0;
	die unless $dat[-1] =~ /^[ATCG]+$/;
	for my $c (@indcnt) {
		my $cnt = 2*$c;
		my $str = substr $dat[-1],0,$cnt;
		my %test;
		++$test{$_} for (split //,$str);
		if ( (sort {$a<=>$b} values %test)[0] == $cnt ) {
			$remained[$c] += $dat[2];
			next;
		}
		my $FH = $OutFH[$c];
		my $value3 = $dat[2] + $remained[$c];
		$remained[$c] = 0;
		print $FH join("\t",@dat[0,1],$value3,$str),"\n";
	}
}

close IN;

close $OutFH[$_] for @indcnt;

__END__
ls -1 20141118_high_coverage.chr*.msmcdat|xargs -n1 ./cutmsmcdat.pl

msmc2 --fixedRecombination -o out.s5.msmc2 chr1s5.msmcdat

msmc --fixedRecombination -o out1/s5 dat/chr*.s5.msmcdat
msmc2 -o out2/s5 dat/chr*.s5.msmcdat

