#!/usr/bin/env perl
use strict;
use warnings;

my ($list,$store,$output) = @ARGV;

my %stat;
open LI,"<$list" or die($!);
open OUT,">$output" or die($!);
while (my $info = <LI>){
	chomp($info);
	my ($M,$F,$C) = split /\s+/,$info;
	my %hete = &get_hete($C);
	my $total = 0;
	my $mismatch = 0;
	open FC,"<$store/p$C.C.tsv" or die($!);
	while (my $line = <FC>){
		chomp($line);
		my @data = split /\t/,$line;
		my @Bases = split /,/,$data[2];
		next unless (defined $hete{$data[0]});
		next if ($data[3] eq '.' or $data[3] < 100);
		my @tmpInfo = splice @data,4;
		my %Dep;
		my $depcheck = 1;
		my $sum = 0;
		for (@tmpInfo){
			my $depsum = 0;
			my @Depinfo = split /[;,]/,$_;
			for my $i (1..scalar @Depinfo - 1){
				if ($Depinfo[$i] eq '.'){$Depinfo[$i] = 0;}
				$Dep{$i - 1} += $Depinfo[$i];
				$depsum += $Depinfo[$i];
			}
			$sum += $depsum;
			if ($depsum < 100){
				$depcheck *= 0;
			}else{
				$depcheck *= 1;
			}
		}
		next if ($depcheck == 0);

		$total++;
                my @dKeys = sort { $Dep{$b} <=> $Dep{$a} } keys %Dep;
                if (@dKeys>1 and $Dep{$dKeys[1]}  >= $Dep{$dKeys[0]} * 0.1){
                        my @rKeys = sort {$a<=>$b} @dKeys[0,1];
                        my $gt = join('/',$Bases[$rKeys[0]],$Bases[$rKeys[1]]);
			unless ($gt eq $hete{$data[0]}){
				$mismatch++;
				print "$C\t$data[0]\t$gt\t$hete{$data[0]}\n";
			}
                }else{
			$mismatch++;
#			print "$C\t$data[0]\n";
		}
	}
	close FC;
	my $freq;
	if ($total == 0){
		$freq = "NA";
	}else{
		$freq = $mismatch/$total;
	}
	print OUT "p$C\t$mismatch/$total\t$freq\n";
}
close LI;

sub get_hete {
	my ($input) = @_;
	my %temp;
	open TE,"<$store/p$input.M.tsv" or die($!);
	while (my $line = <TE>){
		chomp($line);
		my @data = split /\t/,$line;
		next if ($data[3] eq '.' or $data[3] < 100);
		my @Depinfo = split /[;,]/,$data[4];
		my @Bases = split /,/,$data[2];
		my %Dep;
		my $depsum = 0;
		for my $i (1..scalar @Depinfo - 1){
			if ($Depinfo[$i] eq '.'){$Depinfo[$i] = 0;}
			$Dep{$i - 1} = $Depinfo[$i];
			$depsum += $Depinfo[$i];
		}
		next if ($depsum <= 50);
		my @dKeys = sort { $Dep{$b} <=> $Dep{$a} } keys %Dep;
		if (@dKeys>1 and $Dep{$dKeys[1]}  >= $Dep{$dKeys[0]} * 0.1){
			my @rKeys = sort {$a<=>$b} @dKeys[0,1];
			my $gt = join('/',$Bases[$rKeys[0]],$Bases[$rKeys[1]]);
			$temp{$data[0]} = $gt;
		}else{
			next;
		}
	}
	close TE;
	return %temp;
}
