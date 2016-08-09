#!/usr/bin/env perl
use strict;
use warnings;

die "Usage: $0 <SimLst> <AnalyseRes> [Shift=10]\n" if @ARGV <2;

my $total=shift;
my $result=shift;
my $shift=shift;
$shift = 10 unless defined $shift;
open TT,$total or die $!;
open RT,$result or die $!;
my (%hum,%vir,%stat);

while(<TT>){
	chop;
	my @a=split;
	#print $a[3]."\n";
	next unless(/\w/);
	my $tmp = $a[2];
	for ( my $kk=$tmp-$shift;$kk<$tmp+$shift;$kk++) {
		$hum{$kk}=1;
	}
	for my $kk ($a[5] .. $a[6]) {
		++$vir{$kk};
	}

}
close TT;

my $chrid = 'gi|568815597|ref|NC_000001.11|';
while(<RT>){
	chomp;
	my @a=split;
	#print $a[2]."\n";
	my $flag=0;
	if(/\t$chrid\t/){
		for my $i ($a[2] .. $a[3]) {
			$flag |=1 if exists $hum{$i};
		}
		for my $i ($a[6] .. $a[7]) {
			$flag |=2 if exists $vir{$i};
		}
		if($flag){
			print "$flag\t$_\n";
		}
	}
	++$stat{$flag};
}
close RT;

for my $k (sort {$a <=> $b} keys %stat) {
	print "\t$k: $stat{$k}\n";
}

__END__
grep -h \> ../sim90/*.Ref.fa |sed 's/>Ref_/Ref:/g'|sed 's/Vir_/Vir:/'|sed 's/R_/R:/'|sed 's/_/ /g'|cat -n > s90.lst
cat s90_analyse/*.analyse > s90_analyse.lst




./overlap.pl simed.lst simVir4_analyseAll.txt >simgot.lst

./overlap.pl simed.lst simVir4_analyseAll.txt 3|sort|uniq|wc -l
	142
./overlap.pl simed.lst simVir4_analyseAll.txt 4|sort|uniq|wc -l
	178
./overlap.pl simed.lst simVir4_analyseAll.txt 5|sort|uniq|wc -l
	186
./overlap.pl simed.lst simVir4_analyseAll.txt 10|sort|uniq|wc -l
	193
./overlap.pl simed.lst simVir4_analyseAll.txt 20|sort|uniq|wc -l
	200

grep chr18 simVir4_analyseAll.txt|sort|uniq|wc -l
	247
grep RefCut simVir4_analyseAll.txt|sort|uniq|wc -l
	277

grep \> simout_*.Ref.fa | sed 's/^simout_m//'|sed 's/.Ref.fa:>/\t/'|sed 's/Ref_/Ref:/g'|sed 's/Vir_/Vir:/'|sed 's/R_/R:/'|sed 's/_/ /g'|cat -n >simed.lst

