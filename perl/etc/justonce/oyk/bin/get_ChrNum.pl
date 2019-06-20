#!/usr/bin/env perl
use strict;
use warnings;

my ($input,$snp_list,$father) = @ARGV;

my %chr = &get_chr($snp_list);
open IN,"<$input" or die($!);
my $first = <IN>;
my ($cpi,$cpe) = ("NA","NA");
my $cnf = 0;
my (%total,%mis);
while (<IN>){
	chomp;
	my @data =  split /\t/,$_;
	if ($_ =~ /^# CPI:/){
		($cpi) = $_ =~ /# CPI: (\S+)$/;
	}elsif ($_ =~ /^# CPE:/){
		($cpe) = $_ =~ /# CPE: (\S+)$/;
	}elsif ($_ =~ /^#/){
		next;
	}else{
		$total{$chr{$data[0]}}++;
		if ($data[7] == 0.0001){
			$mis{$chr{$data[0]}}++;
			if ($chr{$data[0]} =~ /^\d+$/){
				$cnf++;
			}
		}
	}
}
close IN;
my $result_cpi = sprintf("%.2e", $cpi);

print "CaseNo:\t$father\n";
print "Con:\n";
print "CPI:\t$result_cpi\n";
print "CPE:\t$cpe\n";
print "CNF:\t$cnf\n";
print "CHR\tTotal\tSite\n";
my ($tsum,$ssum) = (0,0);
foreach my $chr (1..22){
	unless (defined $total{$chr}){
		$total{$chr} = 0;
	}
        unless (defined $mis{$chr}){
                $mis{$chr} = 0;
        }
	my $last = $total{$chr} - $mis{$chr};
	$tsum += $total{$chr};
	$ssum += $last;
	print "$chr\t$total{$chr}\t$last\n";
}
print "Sum\t$tsum\t$ssum\n";

###################
sub get_chr {
	my $tempin = shift;
	my %temp;
	open TE,"<$tempin" or die($!);
	while (<TE>){
		chomp;
		my @data = split /\t/,$_;
		$data[1] =~ s/chr//;
		$temp{$data[0]} = $data[1];
	}
	close TE;
	return %temp;
}
