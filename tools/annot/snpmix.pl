#!/bin/env perl
use lib '/share/raid010/resequencing/resequencing/tmp/bin/annotation/glfsqlite';
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.2;

our $opts='i:o:s:f:bv';
our ($opt_i, $opt_o, $opt_s, $opt_v, $opt_b, $opt_f);

our $help=<<EOH;
\t-i PSNP list (./psnp.lst) for chrid.individual.finalSNPs
\t-f fabyChr path (./faByChr) for [chrid].fa
\t-s GLF list (./glf.list), will use \$1 of (([^/]+)/[^/]+$) for sample names
\t-o Output Prefix (./indGenomes/ig_)
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_o='./indGenomes/ig_' if ! defined $opt_o;
$opt_i='./psnp.lst' if ! $opt_i;
$opt_s='./glf.list' if ! $opt_s;
$opt_f='./faByChr' if ! $opt_f;

$opt_i =~ s/\/$//;
$opt_f =~ s/\/$//;
print STDERR "From [$opt_i] to [$opt_o], with [$opt_s][$opt_f]/\n";
if (! $opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];

system('mkdir','-p',$opt_o);
system('rmdir',$opt_o) if $opt_o =~ m#/[\s.]*[^/\s.]+[^/]*$#;

my @Samples;
open L,'<',$opt_s or die "[x]Error opening $opt_s: $!\n";
print STDERR "[!]Sample Order: ";
while (<L>) {
	m#([^/]+)/[^/]+$#;
	push @Samples,$1;
	print STDERR (scalar @Samples),':[',$1,"] ";
}
print STDERR "\n";

print STDERR "[!]Parsing SNP ";
open P,'<',$opt_i or die "[x]Error opening $opt_i: $!\n";
while (my $file=<P>) {
	chomp $file;
	open SNP,'<',$file or (warn "\n[!]Error opening $file: $!\n" and next);
	print STDERR ".\b";
	my (%SNP,$chr,$pos,$ref,$tail,$i);
	while (<SNP>) {
		($chr,$pos,$ref,$tail)=split /\t/;
		my @indSNP;
		for (split / /,$tail) {
			next unless /[ACGTRYMKSWHBVDNX-]/;
			s/-/n/;
			push @indSNP,$_;
		}
		$SNP{$pos}=[$ref,\@indSNP];
	}
	print STDERR "+\b";

	my @FH;
	for (@Samples) {
		$file=$opt_o.$_.'.'.$chr.'.fa';
		my $fh;
		open $fh,'>',$file or die "[x]Error opening $file: $!\n";
		print $fh ">${_}---$chr\n";
		push @FH,$fh;
	}
	warn '[!]PSNP:[',1+$#{${$SNP{$pos}}[1]},'] != File:[',(scalar @FH),"].\n" if $#FH != $#{${$SNP{$pos}}[1]};

	$file=$opt_f.'/'.$chr.'.fa';
	open FA,'<',$file or (warn "\n[!]Error opening $file: $!\n" and next);
	$ref=<FA>;
	$ref =~ />(\S+)/;
	warn "[!]Different ChrID, [>$1] in [$file] !\n" if $1 ne $chr;
	$i=1;
	while ($ref=getc(FA)) {
		unless ($ref =~ /[ACGTRYMKSWHBVDNX]/i) {
			last if $ref eq '>';
			next;
		}
		unless ($i%80) {
			print $_ "\n" for @FH;
		}
		unless ($SNP{$i}) {
			print $_ $ref for @FH;
		} else {
			my ($refbase,$indSNPr)=@{$SNP{$i}};
			warn "[!]RefBase differ, SNP:[$refbase] ne FASTA:[$ref].\n" if $refbase ne uc($ref);
			my $t=0;
			for (@$indSNPr) {
				$tail=$FH[$t];
				print $tail $_;
				++$t;
			}
		}
		++$i;
	}
	print STDERR '-';
}

my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
