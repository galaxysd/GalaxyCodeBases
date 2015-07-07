#!/usr/bin/env perl
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );

use Data::Dump qw(ddx);

die "Usage: $0 <Ref> <Virus> <Outprefix>\n" if @ARGV <3;

my ($Reff,$Virf,$outp)=@ARGV;

$Reff='/share/users/huxs/git/toGit/perl/readsim/Chr1.fa';
$Virf='/share/users/huxs/work/bsvir/HBV.AJ507799.2.fa';

my $SampleCnt = 100;
my $Depth = 50;
my $PEinsertLen=200;
my $VirFragMax = 500;
my $VirFragMin = 20;
my $RefBorder = $PEinsertLen + 1000;
my $RefNratioMax = 0.02;

sub openfile($) {
	my ($filename)=@_;
	my $infile;
	if ($filename=~/.xz$/) {
		open( $infile,"-|","xz -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.gz$/) {
		open( $infile,"-|","gzip -dc $filename") or die "Error opening $filename: $!\n";
	} elsif ($filename=~/.bz2$/) {
		open( $infile,"-|","bzip2 -dc $filename") or die "Error opening $filename: $!\n";
	} else {open( $infile,"<",$filename) or die "Error opening $filename: $!\n";}
	return $infile;
}

sub getRefChr1st($) {
	my $GENOME = $_[0];
	while (<$GENOME>) {
		s/^>//;
		/^(\S+)/ or next;
		my $seqname = $1;
		print STDERR ">$seqname ...";
		$/=">";
		my $genome=<$GENOME>;
		chomp $genome;
		$genome=~s/\s//g;
		$/="\n";
		my $thelength = length $genome;
		print STDERR "\b\b\b", $thelength, ".\n";
		return $genome;
	}
}

my $Refh = openfile($Reff);
my $Refstr = getRefChr1st($Refh);
my $RefLen = length $Refstr;
close $Refh;
my $Virfh = openfile($Virf);
my $Virstr = getRefChr1st($Virfh);
my $VirLen = length $Virstr;
close $Virfh;
$Virstr .= $Virstr;	# circle


my (@Refticks);
while (@Refticks < 100) {
	my $pos0 = int(rand($RefLen-(2*$RefBorder)))+$RefBorder;
	my $str0 = substr $Refstr,($pos0-$PEinsertLen),2*$PEinsertLen;
	my $seq = $str0;
	my $N = $seq=~tr/Nn//;
	next if $N > 2*$PEinsertLen*$RefNratioMax;
	push @Refticks,$pos0
}
@Refticks = sort {$a<=>$b} @Refticks;

ddx \@Refticks;

__END__
./virusinserts.pl /share/users/huxs/git/toGit/perl/readsim/Chr1.fa /share/users/huxs/work/bsvir/HBV.AJ507799.2.fa test
