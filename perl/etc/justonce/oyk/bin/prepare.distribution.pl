#!/usr/bin/env perl
use strict;
use warnings;

my ($list,$store,$output,$depout) = @ARGV;

my (%stat,%dep);
my @fams;
open LI,"<$list" or die($!);
open OUT,">$output" or die($!);
open DOUT,">$depout" or die($!);
while (my $info = <LI>){
	chomp($info);
	next if ($info =~ /^#/);
	my ($M,$F,$C) = split /\s+/,$info;
	my $fam = "p$C";
	push @fams,$fam;
	my %homo = &get_homo($C);
	open FC,"<$store/p$C.C.tsv" or die($!);
	while (my $line = <FC>){
		chomp($line);
		my @data = split /\t/,$line;
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
		$dep{$data[0]}{$fam} = $sum;

		next unless (defined $homo{$data[0]});
		next if ($depcheck == 0);
		my @values = sort {$b<=>$a} values %Dep;
		if (scalar @values > 1){
			my $freq = sprintf "%.1f",$values[1] / $sum * 100;
			$stat{$data[0]}{$fam} = $freq;
#			$stat{$data[0]}{$fam} = int($values[1] / $sum * 100);
		}else{
			$stat{$data[0]}{$fam} = 0;
		}
	}
	close FC;
}
close LI;

print OUT "\t";
print OUT join("\t",@fams),"\n";
foreach my $locus (sort keys %stat){
	print OUT "$locus\t";
	my @values;
	foreach my $fam (@fams){
		if (defined $stat{$locus}{$fam} && $stat{$locus}{$fam} != 0){
			push @values,$stat{$locus}{$fam};
		}else{
			push @values,"";
		}
	}
	print OUT join("\t",@values),"\n";
}
close OUT;

print DOUT "\t";
print DOUT join("\t",@fams),"\n";
foreach my $locus (sort keys %dep){
	print DOUT "$locus\t";
	my @values;
	foreach my $fam (@fams){
		if (defined $dep{$locus}{$fam} && $dep{$locus}{$fam} != 0){
			push @values,$dep{$locus}{$fam};
		}else{
			push @values,"";
		}
	}
	print DOUT join("\t",@values),"\n";
}
close DOUT;

#########################################################
sub get_homo {
	my ($C) = @_;
	my %temp;
	open TE,"<$store/p$C.M.tsv" or die($!);
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

		my @values = sort {$b<=>$a} values %Dep;
		if (scalar @values > 1){
			if ($values[1] <= $values[0] * 0.01){
				$temp{$data[0]}++;
			}
		}else{
			$temp{$data[0]}++;
		}
	}
	close TE;
	return %temp;
}
