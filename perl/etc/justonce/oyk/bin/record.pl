#!/usr/bin/env perl
use strict;
use warnings;

my ($list,$store,$outdir) = @ARGV;

my (%stat,%dep);
my @fams;
open LI,"<$list" or die($!);
while (my $info = <LI>){
	chomp($info);
	next if ($info =~ /^#/);
	my ($M,$F,$C) = split /\s+/,$info;
	&output($C,"M",$M);
	&output($C,"F",$F);
}
close LI;

#########################################################
sub output {
	my ($C,$tag,$id) = @_;
	open TE,"<$store/p$C.$tag.tsv" or die($!);
	open TOUT,">$outdir/$id.txt" or die($!);
	print TOUT "#Sample ID: $id\n";
	print TOUT "MarkerID\tGenoType\tQUAL\tACGTdepths\n";

	my @AllBase = qw{A C G T};
	while (my $line = <TE>){
		chomp($line);
		my @data = split /\t/,$line;
		next if ($data[3] eq '.' or $data[3] < 100);
		my @tM = splice @data,4;
		my %Dep;
		my $depcheck = 1;
		my @bases = split /,/,$data[2];
		for (@tM){
			my $depsum = 0;
			my @Depinfo = split /[;,]/,$_;
			for my $i (1..scalar @Depinfo - 1){
				if ($Depinfo[$i] eq '.'){$Depinfo[$i] = 0;}
				$Dep{$bases[$i - 1]} += $Depinfo[$i];
				$depsum += $Depinfo[$i];
			}
			if ($depsum <= 50){
				$depcheck *= 0;
			}else{
				$depcheck *= 1;
			}
		}
		next if ($depcheck == 0);

		my @sort = sort {$Dep{$b}<=>$Dep{$a}} keys %Dep;
		my $genotype;
		if (scalar @sort == 1){
			$genotype = join "/",$sort[0],$sort[0];
		}elsif (scalar @sort > 1 && $Dep{$sort[1]} <= $Dep{$sort[0]} * 0.01){
			$genotype = join "/",$sort[0],$sort[0];
		}elsif (scalar @sort > 1 && $Dep{$sort[1]} >= $Dep{$sort[0]} * 0.1){
			my @rKeys = sort @sort[0,1];
			$genotype = join "/",@rKeys;
		}

		for (@AllBase){
			unless (defined $Dep{$_}){
				$Dep{$_} = ".";
			}
		}
		if (defined $genotype && $genotype =~ /[ACGT]+/){
			print TOUT "$data[0]\t$genotype\t$data[3]\t$Dep{A},$Dep{C},$Dep{G},$Dep{T}\n";
		}
	}
	close TE;
	close TOUT;
}
