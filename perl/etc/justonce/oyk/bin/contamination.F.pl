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

my $list = shift;
my $store = shift;
my $output = shift;

my @fams;
my %Fs;
open LI,"<$list" or die($!);
open OUT,">$output" or die($!);
while (my $info = <LI>){
	chomp($info);
	my ($M,$F,$C) = split /\s+/,$info;
	push @fams,"p$C";
	$Fs{$M} = "p$C.M";
	$Fs{$F} = "p$C.F";
	$Fs{$C} = "p$C.C";
}
close LI;

my @fathers = sort keys %Fs;
print OUT "\t";
print OUT join("\t",@fathers),"\n";
foreach my $family (@fams){
	my $Mfile = "$store/$family.M.tsv";
	my $Cfile = "$store/$family.C.tsv";
	print OUT "$family\t";
	my @outputs;
	foreach my $father (@fathers){
		my $total = 0;
		my $mismatch = 0;
		my $Ffile = "$store/$Fs{$father}.tsv";
		system("mkdir -p $store/temp") unless (-e "$store/temp");
		print "perl $RealBin/oykn.pl $thetMode $thetParentage $Mfile $Ffile $Cfile $store/temp/tm$family\n";
		my @info = readpipe("perl $RealBin/oykn.pl $thetMode $thetParentage $Mfile $Ffile $Cfile $store/temp/tm$family");
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
		my $outinfo = join "/",$mismatch,$total;
		push @outputs,$outinfo;
	}
	print OUT join("\t",@outputs),"\n";
}
close OUT;

system("rm -rf $store/temp") if (-e "$store/temp");
