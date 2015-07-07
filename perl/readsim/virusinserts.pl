#!/usr/bin/env perl
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );

use Data::Dump qw(ddx);

die "Usage: $0 <Ref> <Virus> <Outprefix>\n" if @ARGV <3;

my ($Reff,$Virf,$outp)=@ARGV;

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
	my $genome;
	while (<$GENOME>) {
		s/^>//;
		/^(\S+)/ or next;
		my $seqname = $1;
		print STDERR " >$seqname ...";
		$/=">";
		my $genome=<$GENOME>;
		chomp $genome;
		$genome=~s/\s//g;
		$/="\n";
		my $thelength = length $genome;
		print STDERR "\b\b\b", $thelength, ".\n";
	}
	return $genome;
}

my $Refh = openfile($Reff);
my $Refstr = getRefChr1st($Refh);
close $Refh;
my $Virfh = openfile($Virf);
my $Virstr = getRefChr1st($Virfh);
close $Virfh;



__END__
./virusinserts.pl /share/users/huxs/git/toGit/perl/readsim/Chr1.fa /share/users/huxs/work/bsvir/HBV.AJ507799.2.fa test
