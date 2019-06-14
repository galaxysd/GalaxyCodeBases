#!/usr/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 1.0.0 @ 20180614
=cut
use strict;
use warnings;
use Cwd qw(abs_path cwd);
use Parse::CSV;
#use lib '.';
use Data::Dump qw(ddx);

my ($bedfn,$cgLfn,$cgRfn,$outfn) = @ARGV;
if (@ARGV < 4) {
	die "Usage: $0 <bed> <Left.cg> <Right.cg> <out.tsv>\n";;
}

my $cgLin = Parse::CSV->new(file => $cgLfn, names => 1, sep_char => "\t");
my $cgRin = Parse::CSV->new(file => $cgRfn, names => 1, sep_char => "\t");

open B,'<',$bedfn or die $?;
my (%Dbed,%Rbed);
while (<B>) {
	my @L = split /\t/;
	next unless $L[2];
	push @{$Dbed{$L[0]}},[$L[1],$L[2]];
	my @rng = $L[1] .. $L[2];
	for (@rng) {
		++$Rbed{$L[0]}{$_};
	}
}
#ddx %Dbed;

sub loadCG($$) {
	my ($rngD,$csvh) = @_;
	my %Dcg;
	while ( my $value = $csvh->fetch ) {
		#ddx $value;
		if (exists $rngD->{ $value->{'#CHROM'} }{ $value->{'POS'} }) {
			my ($cgCnt,$cgC,$cgW)=(0,0,0);
			if ($value->{'Crick-COVERAGE'} ne '.') {
				if ($value->{'Crick-COVERAGE'} >= 5) {
					$cgC = $value->{'Crick-METH'}/$value->{'Crick-COVERAGE'};
					$cgCnt += 1;
				}
			}
			if ($value->{'Watson-COVERAGE'} ne '.') {
				if ($value->{'Watson-COVERAGE'} >= 5) {
					$cgW = $value->{'Watson-METH'}/$value->{'Watson-COVERAGE'};
					$cgCnt += 1;
				}
			}
			if ($cgCnt) {
				$Dcg{$value->{'#CHROM'}}{$value->{'POS'}} = [$cgCnt/2, ($cgC+$cgW)/$cgCnt];
			}
		} else {
			next;
		}
	}
	return \%Dcg;
}

my %DcgL = %{ loadCG(\%Rbed,$cgLin) };
ddx \%DcgL;
my %DcgR = %{ loadCG(\%Rbed,$cgRin) };
ddx \%DcgR;



__END__
./cdmr.pl gencode.v30.annotation.bed.h B7B.cg.h Normal.cg.h o.tsv
