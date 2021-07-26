#!/usr/bin/env perl

use strict;
use warnings;
#use Data::Dump qw(ddx);

my $vcf = shift;
my $snpLst = shift;
my $minQual = shift;
my $minDepth = shift;

$minQual = 0 unless (defined $minQual);
$minDepth = 100 unless (defined $minDepth);

my (%need,%record);
&read_list(\%need,$snpLst);
open IN,"<$vcf" or die($!);
while (my $line = <IN>){
        chomp($line);
        next if ($line =~ /^#/);
        my @data = split /\t/,$line;
        next unless (defined $need{$data[2]});
	next if ($data[5] < $minQual);
        my $depth;
	if ($data[7] =~ /DP=([0-9]+);/){
		$depth = $1;
	}else{
		$depth = 0;
	}

        my %pos;
        my @info = split /:/,$data[8];
        for my $i (0..$#info){
                $pos{$info[$i]} = $i;
        }
        if ($depth < $minDepth){
	        $record{$data[2]} = "$depth\t$data[5]\tUnknown";
	}else{
                my $bases = join ",",$data[3],$data[4];
                my @allele = split /,/,$bases;
                my @infov = split /:/,$data[9];
                my @AD = split /,/,$infov[$pos{'AD'}];
                my %temp;
                for my $i (0..$#AD){
	                $temp{$allele[$i]} = $AD[$i];
                }
                my @sortA = sort {$temp{$b} <=> $temp{$a}} keys %temp;
                if (!defined $sortA[1]){
                        $record{$data[2]} = join "\t",$depth,$data[5],"$sortA[0]/$sortA[0]";
                }elsif($temp{$sortA[1]} / $temp{$sortA[0]} >= 0.1){
                        $record{$data[2]} = join "\t",$depth,$data[5],"$sortA[0]/$sortA[1]";;
                }elsif ($temp{$sortA[1]} / $temp{$sortA[0]} < 0.02){
                        $record{$data[2]} = join "\t",$depth,$data[5],"$sortA[0]/$sortA[0]";
                }else{
	                $record{$data[2]} = "$depth\t$data[5]\tUnknown";;
	        }
	}       
}       
close IN;

foreach my $snp (sort {$need{$a}->[0] <=> $need{$b}->[0]} keys %need){
	if (defined $record{$snp}){
		print "$snp\t$record{$snp}\t";
		print join(":",$need{$snp}->[1],$need{$snp}->[2]),"\n";
	}else{
		print "$snp\t0\t0\tUnknown\t";
		print join(":",$need{$snp}->[1],$need{$snp}->[2]),"\n";
	}
}

###############
sub read_list {
        my ($hash,$tempIN) = @_;
	my $te;
	if ($tempIN =~ /gz$/){
		open $te,"zcat $tempIN|" or die($!);
	}else{
        	open $te,"<$tempIN" or die($!);
	}
        while (<$te>){
                chomp;
		next if ($_ =~ /^#/);
                my @data = split /\t/,$_;
		$hash->{$data[3]} = [$.,$data[0],$data[1]];
	}
	close $te;
}
