#!/usr/bin/env perl
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );

use Data::Dump qw(ddx);

my $SAMTOOLSBIN="samtools";

die "Usage: $0 <output> <sam/bam files>\n" if @ARGV <2;

my $out = shift;

warn "[!]Stat [$out] <- [",join('],[',@ARGV),"]\n";
open OUT,'>',$out or die "Error opening $out: $!\n";

my %StatDat;
for my $filename (@ARGV) {
	print STDERR "[!]Reading [$filename]: ";
	open INSAM,'-|',"$SAMTOOLSBIN view $filename" or die "Error opening $filename: $!\n";
	while (<INSAM>) {
		my @read1=split /\t/;
		my $OPT = join("\t",@read1[11 .. $#read1]);
		# next if $OPT =~ /\b(XT:A:R\b)|(XA:Z:)/;	# bwa mem use XA:Z for Alternative hits, SA:Z for Chimeric reads.
		my @cigar = $read1[5] =~ /(\d+)(\w)/g;
		my ($reflen,$maxM,$readlen)=(0,0,0);
		while (@cigar) {
			my ($len,$op) = splice(@cigar,0,2);
			if ($op eq 'M') {
				$reflen += $len;
				$readlen += $len;
				$maxM = $len if $maxM < $len;
				next;
			}
			$reflen += $len if $op eq 'D';
			$readlen += $len if $op eq 'I';
		}
		$StatDat{$filename}{'Bases'} += $readlen;
		++$StatDat{$filename}{'Reads'};
	}
	close INSAM;
	warn "done.\n";
}

for my $filename (sort keys %StatDat) {
	$StatDat{'/__All__/'}{$_} += $StatDat{$filename}{$_} for keys %{$StatDat{$filename}};
}
my @StatKeys = sort keys %{$StatDat{'/__All__/'}};

my $tmpline = join("\t",'FileName',@StatKeys);
warn "$tmpline\n";
print OUT "#$tmpline\n";
for my $filename (sort keys %StatDat) {
	print STDERR $filename;
	print OUT $filename;
	for my $k (@StatKeys) {
		print STDERR "\t",$StatDat{$filename}{$k};
		print OUT "\t",$StatDat{$filename}{$k};
	}
	warn "\n";
	print OUT "\n";
}
#ddx \%StatDat;

__END__
./bam.pl stat.log  /share/users/huxs/work/tiger2/part1/PD2M-LSJ-BHX01*
