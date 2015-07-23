#!/usr/bin/env perl
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );

use Data::Dump qw(ddx);

die "Usage: $0 <Host> <Virus> <bam sort-n>\n" if @ARGV <3;

my ($Reff,$Virf,$bamin)=@ARGV;

$Reff='/share/users/huxs/git/toGit/perl/readsim/Chr1.fa';
$Virf='/share/users/huxs/work/bsvir/HBV.AJ507799.2.fa';

sub openfile($) {
	my ($filename)=@_;
	my $infile;
	if ($filename=~/.[bs]am(.gz)?$/) {
		open( $infile,"-|","samtools view -F 256 $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.gz$/) {
		open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.bz2$/) {
		open( $infile,"-|","bzip2 -dc $filename") or die "Error opening $filename: $!\n";
	} else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
	return $infile;
}

sub getRefChr1stID($) {
	my $GENOME = $_[0];
	while (<$GENOME>) {
		s/^>//;
		/^(\S+)/ or next;
		my $seqname = $1;
		print STDERR ">$seqname ...";
		return $seqname;
	}
}

my $Refh = openfile($Reff);
my $RefID = getRefChr1stID($Refh);
close $Refh;
my $Virfh = openfile($Virf);
my $VirID = getRefChr1stID($Virfh);
close $Virfh;

my $bamfh = openfile($bamin);
while (<$bamfh>) {
	my @dat1 = split /\t/;
	$_=<$bamfh>;
	my @dat2 = split /\t/;
	die "[x]Read1 & Read2 not match ! [$dat1[0]] ne [$dat2[0]]\n" if $dat1[0] ne $dat2[0];
	warn "[$dat1[0]] [$dat2[0]] [$dat1[3]] [$dat2[3]]\n";
}
close $bamfh;

__END__
./bamstat.pl /share/users/huxs/git/toGit/perl/readsim/Chr1.fa /share/users/huxs/work/bsvir/HBV.AJ507799.2.fa /share/users/huxs/work/bsvir/bsI/Test1_aln/T.sn.bam
