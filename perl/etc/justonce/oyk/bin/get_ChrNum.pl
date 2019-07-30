#!/usr/bin/env perl
use strict;
use warnings;

my ($input,$snp_list,$father) = @ARGV;

my (%chr,%changeID);
&get_chr(\%chr,\%changeID,$snp_list);
open IN,"<$input" or die($!);
my $first = <IN>;
my ($cpi,$cpe) = ("NA","NA");
my $cnf = 0;
my (%total,%mis);
my @MergeoOut;
push @MergeoOut,"检测位点编号\t染色体\t父本基因型\t母本基因型\t胎儿基因型\t是否错配";
while (<IN>){
	chomp;
	my @data =  split /\t/,$_;
	for (@data){
		if ($_ =~ /;\d+,/){
			$_ =~ s/;(\S+)//;
		}
	}
	if ($_ =~ /^# CPI:/){
		($cpi) = $_ =~ /# CPI: (\S+)$/;
	}elsif ($_ =~ /^# CPE:/){
		($cpe) = $_ =~ /# CPE: (\S+)$/;
	}elsif ($_ =~ /^#/){
		next;
	}else{
		$total{$chr{$data[0]}}++;
		push @MergeoOut,"$changeID{$data[0]}\tchr$chr{$data[0]}\t$data[5]\t$data[4]\t$data[6]\t匹配";
		if ($data[7] == 0.0001){
			push @MergeoOut,"$changeID{$data[0]}\tchr$chr{$data[0]}\t$data[5]\t$data[4]\t$data[6]\t错配";
			$mis{$chr{$data[0]}}++;
			if ($chr{$data[0]} =~ /^\d+$/){
				$cnf++;
			}
		}
	}
}
close IN;
my $result_cpi = sprintf("%.2e", $cpi);

my @cpi_info = split /[Ee]/,$result_cpi;
my @cpe_info = split /[Ee]/,$cpe;
my $con;
if ($cpi_info[scalar @cpi_info - 1] >= 4 && $cpe_info[scalar @cpe_info - 1] <= -8){
	$con = "Y";
}elsif ($cpi_info[scalar @cpi_info - 1] <= -4 && $cpe_info[scalar @cpe_info - 1] <= -8){
	$con = "N";
}else{
	die("#Not up to standard.\n");
}

print "CaseNo:\t$father\n";
print "Con:\t$con\n";
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
	if ($con eq 'Y'){
		print "$chr\t$total{$chr}\t$last\n";
	}else{
		print "$chr\t$total{$chr}\t",$total{$chr} - $last,"\n";
	}
}

if ($con eq 'Y'){
	print "Sum\t$tsum\t$ssum\n\n";
}else{
	print "Sum\t$tsum\t",$tsum - $ssum,"\n\n";
}

print join("\n",@MergeoOut),"\n";

###################
sub get_chr {
	my ($hash1,$hash2,$tempin) = @_;
	my $count = 1;
	open TE,"<$tempin" or die($!);
	while (<TE>){
		chomp;
		my @data = split /\t/,$_;
		$data[1] =~ s/chr//;
		$hash1->{$data[0]} = $data[1];
		my $len = 4 - length($count);
		my $tag = join "",0 x $len,$count;
		$hash2->{$data[0]} = "SNP$tag";
		$count++;
	}
	close TE;
}
