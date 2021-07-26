#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use Cwd 'abs_path';

my ($template,$item)= @ARGV;

my %fill;
&get_item(\%fill,$item);

local $/ = ">INSERT<";
open TE,"<$template" or die($!);
while (<TE>){
	$_ =~ s/>INSERT<$//;
	if ($_ eq 'snp table'){
		my $head = "| tag | rsID | pos | genotype |  \n| :---- | :------- | :------- | :---- |";
		my $SNPtxt = $fill{'snp table'};
		my @trans = &mtable($SNPtxt);
		unshift @trans,$head;
		print join("  \n",@trans);
	}elsif ($_ eq 'qual plot'){
		my $qualPlot = $fill{'qual plot'};
		my @file = split /,/,$qualPlot;
		foreach my $plot (@file){
			my $qPath = abs_path($plot);
			my $decode = `base64 $qPath`;
			chomp($decode);
			$decode =~ s/\n//g;
			print "![png qual](data:image/png;base64,$decode)";
		}
	}elsif ($_ eq 'quality control'){
		my $qcData = $fill{'quality control'};
		my @file = split /,/,$qcData;
		my $head = "| lib | bases | gc | low quality |  \n| :---- | :------- | :------- | :---- |";
		my @Qtable;
		push @Qtable,$head;
		foreach my $qcfile (@file){
			my $qc = &read_qc($qcfile);
			push @Qtable,"| $qc->[0] | $qc->[1] | $qc->[2] | $qc->[3] |";
		}
		print join("  \n",@Qtable);
	}elsif ($_ eq 'sample id'){
		print "$fill{'sample id'}";
	}elsif ($_ eq 'library number'){
		print "$fill{'library number'}";
	}elsif ($_ eq "total reads"){
		my $statsF = $fill{'align stats'};
		my $stats = &bwa_stats($statsF);
		print "$stats->[0]";
	}elsif ($_ eq "mapping rate"){
		my $statsF = $fill{'align stats'};
		my $stats = &bwa_stats($statsF);
		print "$stats->[1]";
	}else{
		print "$_";
	}
}
close TE;

########################################
sub get_item {
	my ($hash,$temp) = @_;
	open TE,"<$temp" or die($!);
	while (<TE>){
		chomp;
		my @data = split /\t/,$_;
		$hash->{$data[0]} = $data[1];
	}
	close TE;
}

sub bwa_stats {
	my $temp = shift;
	unless (-e $temp){
		return['NA','NA'];
	}
	local $/ = "\n";
	my ($total,$map) = (1,0);
	open FE,"<$temp" or die($!);
	while (<FE>){
		chomp;
		next unless ($_ =~ /^SN/);
		my @data = split /\t/,$_;
		if ($data[1] eq 'raw total sequences:'){
			$total = $data[2];
		}elsif ($data[1] eq 'reads mapped:'){
			$map = $data[2];
		}
	}
	close FE;
	my $mrate = join " ",sprintf("%.2f",$map/$total*100),'%';
	return[$total,$mrate];
}

sub mtable {
	my $temp = shift;
	unless (-e $temp){
		return ("|    |    |    |");
	}
	my @code = ();
	local $/ = "\n";
	open FE,"<$temp" or die($!);
	while (<FE>){
		chomp;
		my @data = split /\t/,$_;
		if ($data[4] =~ /^chr[XYM]/){
			push @code,"| sexual | $data[0] | $data[4] | $data[3] |";
		}else{
			push @code,"| identity | $data[0] | $data[4] | $data[3] |";
		}
	}
	close FE;
	return @code;
}

sub read_qc {
	my $temp = shift;
	unless (-e $temp){
		return ("|    |    |    |");
	}
	my ($libnm,$rnum,$len,$gcp,$low);
	local $/ = "\n";
	open FE,"<$temp" or die($!);
	while (<FE>){
		chomp;
		my @data = split /\t/,$_;
		if ($data[0] eq 'Filename'){
			$libnm = $data[1];
			$libnm =~ s/\.fq\.gz$//;
		}elsif ($data[0] eq 'Total Sequences'){
			$rnum = $data[1];
		}elsif ($data[0] eq 'Sequences flagged as poor quality'){
			$low = $data[1];
		}elsif ($data[0] eq 'Sequence length'){
			$len = $data[1];
		}elsif ($data[0] eq '%GC'){
			$gcp = "$data[1]%";
		}
	}
	close FE;
	my $bases = join " ",sprintf("%.2f", $rnum * $len / 1000000),'M';
	return [$libnm,$bases,$gcp,$low];
}
