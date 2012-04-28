#!/bin/env perl
use lib 'E:/BGI/toGit/perlib/etc';
use strict;
use warnings;
use File::Basename;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.1;


our $opts='p:b';
our($opt_p, $opt_b);

our $help=<<EOH;
\t-p output prefix
\t-b No pause for batch runs
EOH

ShowHelp();

die "[x]Must specify output prefix !\n" if ! $opt_p;

my $FASTAwidth=75;

open O,'>',$opt_p.'.fa' or die "Error: $!\n";

sub openfile($) {
	my $opti =$_[0];
	my $infile;
	if ($opti =~ /.bz2$/) {
		open( $infile,"-|","bzip2 -dc $opti") or die "Error: $!\n";
	} elsif ($opti =~ /.gz$/) {
		open( $infile,"-|","gzip -dc $opti") or die "Error: $!\n";
	} else {open( $infile,"<",$opti) or die "Error: $!\n";}
	return $infile;
}

warn "Files: [@ARGV]\n";

my $TotalLen=0;
while ($_=shift @ARGV) {
	warn "[$_]\n";
	my $infile=openfile($_);
	{
		local $/=">";
		$_=<$infile>;
		die "[x]Not a FASTA file !\n" unless /^\s*>/;
	}
	while (<$infile>) {
		chomp;
		my $Head;
		my ($id,$desc)=split / /,$_,2;
		if ($desc && $desc !~ /^\s*$/) {
			$desc=~s/\t/_/g;
			$Head="$id $desc";
		} else { $desc='.';$Head=$id; }
		$/=">";
		my $seq=<$infile>;
		chomp $seq;
		$seq =~ s/\s//g;
		$/="\n";
		my $len=length($seq);
		my $lc = $seq=~tr/agct/nnnn/;
		$TotalLen += $len;
		print STDERR ">$id\t$len, $lc\t[$desc]:\n";
		print O ">$Head $desc\n";
		for (my $i=0; $i<$len; $i+=$FASTAwidth) {
			print O substr($seq,$i,$FASTAwidth)."\n";
		}
	}
	close $infile;
}
close O;

my $cmd='bwa index ';
if ($TotalLen > 2000000000) {
	$cmd .= '-a bwtsw ';
}
system($cmd."-p $opt_p ${opt_p}.fa");
system("samtools faidx ${opt_p}.fa");
