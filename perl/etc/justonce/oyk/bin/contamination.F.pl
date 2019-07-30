#!/usr/bin/env perl
=pod
Author: HU Xuesong @ BGI <huxuesong@genomics.org.cn>, LI Bowen <libowen@genomics.cn>
Version: 1.0.0 @ 20180624
=cut
use strict;
use warnings;
use POSIX;

use FindBin qw($RealBin);
if ($FindBin::VERSION < 1.51) {
	warn "[!]Your Perl is too old, thus there can only be ONE `bsuit` file in your PATH. [FindBin Version: $FindBin::VERSION < 1.51]\n\n"
}
FindBin::again();
use lib "$RealBin/../";

my @tModes = qw(CHIP PCR);
my %tMode = map { $_ => 1 } @tModes;
my @tParentages = qw(DUO TRIO);
my %tParentage = map { $_ => 1 } @tParentages;

my $thetMode = uc shift;
unless (exists $tMode{$thetMode}) {
	die "[x]mode can only be:[",join(',',@tModes),"].\n";
}
my $thetParentage = uc shift;
unless (exists $tParentage{$thetParentage}) {
	die "[x]Parentage can only be:[",join(',',@tParentages),"].\n";
}

my $qc_dir = shift;
my $output = shift;
open OUT,">$output" or die($!);

$qc_dir =~ s/\/$//g;
my @path = split /\//,$qc_dir;
my $now = $path[scalar @path - 2];
my @useless = splice @path,scalar @path - 2;
my $total_dir = join "/",@path;
my @content = readpipe("ls -F $total_dir");

my %chips;
for (@content){
	chomp($_);
	if ($_ =~ /^\d+CG\S+\/$/ or $_ =~ /^\d+CG\S+\@$/){
		$_ =~ s/\@$//;
		$_ =~ s/\/$//;
		my ($date) = $_ =~ /^(\d+)CG/;
		$chips{$_} = $date;
	}
}

my @dirs = ($now);
if ($now =~ /^\d+CG/){
	my ($nowtime) = $now =~ /^(\d+)CG/;
	foreach my $dir (sort {$chips{$b}<=>$chips{$a}} keys %chips){
		next unless (-e "$total_dir/$dir/6record");
		next if ($dir eq $now);
		if ($chips{$dir} eq $nowtime){
			push @dirs,$dir;
		}elsif ($chips{$dir} < $nowtime && scalar @dirs < 3){
			push @dirs,$dir;
		}
	}
}

my @fams;
my (%Fs,%store);
foreach my $take_dir (@dirs){
	my $list = "$total_dir/$take_dir/family.lst";
	open LI,"<$list" or die($!);
	while (my $info = <LI>){
		chomp($info);
		my ($M,$F,$C) = split /\s+/,$info;
		if ($take_dir eq $now){
			push @fams,"p$C";
		}
		$Fs{$M} = "p$C.M";
		$Fs{$F} = "p$C.F";
		$store{$M} = "$total_dir/$take_dir/4tsv";
		$store{$F} = "$total_dir/$take_dir/4tsv";
	}
	close LI;
}

print join("\n",@fams),"\n";
my @fathers = sort keys %Fs;
print OUT "\tTotal\t";
print OUT join("\t",@fathers),"\n";
foreach my $family (@fams){
	my $Mfile = "$total_dir/$now/4tsv/$family.M.tsv";
	my $Cfile = "$total_dir/$now/4tsv/$family.C.tsv";
	print OUT "$family\t";
	my @outputs;
	my @totals;
	foreach my $father (@fathers){
		my $total = 0;
		my $mismatch = 0;
		my $Ffile = "$store{$father}/$Fs{$father}.tsv";
		system("mkdir -p $total_dir/$now/4tsv/temp") unless (-e "$total_dir/$now/4tsv/temp");
		my @info = readpipe("perl $RealBin/oykn.pl $thetMode $thetParentage $Mfile $Ffile $Cfile $total_dir/$now/4tsv/temp/tm$family.$father");
		foreach my $line (@info){
			chomp($line);
			next if ($line =~ /^#/);
			my @data = split /\t/,$line;
			next unless (defined $data[7]);
			$total++;
			if ($data[7] == 0.0001){
				$mismatch++;
			}
		}
		push @outputs,$mismatch;
		push @totals,$total;
	}
	my @sort = sort {$a<=>$b} @totals;
	print OUT "$sort[0]-$sort[scalar @sort - 1]\t";
	print OUT join("\t",@outputs),"\n";
}
close OUT;

system("rm -rf $total_dir/$now/4tsv/temp") if (-e "$total_dir/$now/4tsv/temp");
