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
my %hash;

while(<TT>){
	chop;
	my @a=split;
	#print $a[3]."\n";
	next unless(/\w/);
	for( my $kk=$a[3]-$shift;$kk<$a[3]+$shift;$kk++){
		$hash{$kk}=1;
	}

}
close TT;

sub docheck($$$) {
	my ($chr,$start,$end) = @_;
	my $ret = 0;
	($start,$end) = sort {$a <=> $b} ($start,$end);
	for my $p ($start .. $end) {
		if (exists $hash{$p}) {
			$ret = 1;
			last;
		}
	}
	return $ret;
}

while(<RT>){
	chomp;
	my @a=split;
	my ($chr,$start,$end)=split /_/,$a[3];
	my $ret = docheck($chr,$start,$end);
	print "$_\n" if $ret;
}
close RT;

__END__
grep \> simout_*.Ref.fa | sed 's/^simout_m//'|sed 's/.Ref.fa:>/\t/'|sed 's/Ref_/Ref:/g'|sed 's/Vir_/Vir:/'|sed 's/R_/R:/'|sed 's/_/ /g'|cat -n >simed.lst

$ wc VirusFinder2/results-virus-loci.txt
     630    9450   58915 VirusFinder2/results-virus-loci.txt

$ ./overlapVF.pl simed.lst VirusFinder2/results-virus-loci.txt |wc
     408    6120   38407
