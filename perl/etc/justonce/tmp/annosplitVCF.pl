#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw (ddx);

my $gtff = '/bak/seqdata/genomes/TigerRefGenome/P_tigris.gene.gtf';
my $wcl = 374233;
#$gtff = 't.gtf';

my ($lastStrand,%GTF,%tcds,$lastchr,%OUTF);

open I,'-|',"gzip -dc $ARGV[0]" or die;
my @Types = ('I',1,2,3,'N');
for (@Types) {
	open $OUTF{$_},'|-',"gzip -9c >vcfout.$_.vcf.gz";
}

while (<I>) {
	if (/^#/) {
		for my $i (@Types) {
			my $fh = $OUTF{$i};
			print $fh $_;
		}
		if (/^##contig=<ID=(\w+),length=(\d+)>$/) {
			$GTF{$1}=' ' x $2;
			print ">$1: $2\n";
		}
		last if /^#CHROM/;
		next;
	}
	die "Err $_";
}

open G,'<',$gtff;
my $tl = 0;
while (<G>) {
	my ($chr,undef,$type,$start,$end,undef,$strand) = split /\t/;
	++$tl;
	next unless exists $GTF{$chr};
	if ($type eq 'transcript') {
		my $str = $GTF{$chr};
		my $len = $end-$start;
		substr $GTF{$chr},$start-1,$len,'I'x$len;

		if (defined $lastStrand) {
			my @Poses = keys %{$tcds{$chr}};
			if ($lastStrand eq '+') {
				@Poses = sort { $a <=> $b } @Poses;
			} else {
				@Poses = sort { $b <=> $a } @Poses;
			}
			my $step = 1;
			for (@Poses) {
				substr $GTF{$chr},$_-1,1,"$step";
				++$step;
				$step = 1 if $step > 3;
			}
			%tcds = ();
		}
		$lastStrand = undef;
		$lastchr = undef;
	} elsif ($type eq 'exon') {
		for ($start .. $end) {
			$tcds{$chr}{$_} = 'E';
		}
		$lastStrand = $strand;
		$lastchr = $chr;
	}
	print STDERR int(100*$tl/$wcl)/100,"%      $tl      \r" unless $tl % 1000;
}
if (defined $lastStrand) {
	my @Poses = keys %{$tcds{$lastchr}};
	if ($lastStrand eq '+') {
		@Poses = sort { $a <=> $b } @Poses;
	} else {
		@Poses = sort { $b <=> $a } @Poses;
	}
	my $step = 1;
	for (@Poses) {
		substr $GTF{$lastchr},$_-1,1,"$step";
		++$step;
		$step = 1 if $step > 3;
	}
	%tcds = ();
}
warn "Done !      \n";
close G;

#print unpack('c',$GTF{'scaffold53'}{44227}),"|\n";
#print $GTF{'scaffold53'}{44227},"|\n";
#ddx \%GTF;die;
open X,'>','gtf.db';
for (sort keys %GTF) {
	print X join("\t",$_,$GTF{$_}),"\n";
}
close X;
warn "OK !      \n";

while (<I>) {
	my ($chr,$pos) = split /\t/;
	my $str = $GTF{$chr};
	my $type = substr $str,$pos-1,1;
	if ($type ne ' ') {
		my $fh = $OUTF{$type};
		print $fh $_;
	} else {
		my $fh = $OUTF{'N'};
		print $fh $_;
	}
}

ddx \%GTF;
for (@Types) {
	close $OUTF{$_};
}

