#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dump qw(ddx);

# http://www.perlmonks.org/?node_id=393426
use FindBin;
use lib $FindBin::Bin;
use getfuse;

our $DEBUG=2.9;
our $minLen = 5;

die "Usage: $0 <ins_size> <hum_sorted_bam> <vir_sorted_bam> [out_prefix]\n" if @ARGV < 3;
my ($isize,$inhum,$invir,$out)=@ARGV;

unless (defined $out) {
	$out='fuse_'.$inhum;
	$out =~ s/\.[sb]am(\.gz)?//g;
}
warn "[!]IN: h=[$inhum],v=[$invir],i=[$isize] to [$out].*\n";

our (%Genome,%ChrLen);	# Let's share them.

my ($ret,$finHum,$finHBV);
$finHum = getsamChrLen($inhum);
$finHBV = getsamChrLen($invir);

#ddx \%ChrLen;
ddx $finHum;
ddx $finHBV;

warn "[!]Ref: h=[$$finHum[0]],v=[$$finHBV[0]]\n";

my $SearchPanRange = 2*$isize + 100;
sub doPileUp() {
	my ($t)=@_;
	my (%PileUp,%Depth);
}
sub doPan() {
	my ($t)=@_;
	my (%PileUp,%Depth);
}


open INVIR,"-|","samtools view $invir" or die "Error opening $invir: $!\n";
my (%PileUpVIR,%DepthVIR);
while (<INVIR>) {
	chomp;
	my ($id, $flag, $ref, $pos, $mapq, $CIAGR, $mref, $mpos, $isize, $seq, $qual, @OPT) = split /\t/;
	my $mappedSeq = parseCIGAR($seq,$CIAGR);
	my $mappedQual = parseCIGAR($qual,$CIAGR);
	my $mappedLen = length $mappedSeq;
	if ( $DEBUG > 1 ) {
		push @{$PileUpVIR{$ref}},[$pos,$mappedLen,$id,$mappedSeq,$mappedQual,$seq,$CIAGR];
	} else {
		push @{$PileUpVIR{$ref}},[$pos,$mappedLen,$id,$mappedSeq,$mappedQual];
	}
}
for my $chr (keys %PileUpVIR) {
	my $thischrlen = $ChrLen{$chr};
	warn "[!]Vir: $chr $thischrlen:\n";
	$DepthVIR{$chr}->[$thischrlen]=0;
	for my $Datp ( @{$PileUpVIR{$chr}} ) {
		for my $pos ( ($Datp->[0]) .. ($Datp->[0] + $Datp->[1]) ) {
			++$DepthVIR{$chr}->[$pos];
		}
	}
	for (@{$DepthVIR{$chr}}) {
		$_ = '.' if not defined $_;
	}
	#warn "@{$DepthVIR{$chr}}\n";
}


__END__
./getfuse.pl 200 grep_s00_C.bs.virsam.sort.bam grep_s00_C.bs.sort.bam
记得把流程中的文件名改下，突出ref是哪个物种。


CMD Log:
bsmap -u -d HomoGRCh38.fa -a F12HPCCCSZ0010_Upload/s00_C.bs_1.fq.gz -b F12HPCCCSZ0010_Upload/s00_C.bs_2.fq.gz -z 64 -p 12 -v 10 -q 2 -o s00_C.bs.bam >s00_C.bs.log 2>s00_C.bs.err

./getUnPaired.pl s00_C.bs.bam	# 出UnPaired的fq和sam，及grep_s00_C.insertsize。

bsmap -z 64 -p 12 -v 10 -q 2 -d HBV.AJ507799.2.fa -a grep_s00_C.bs.1.fq.gz -b grep_s00_C.bs.2.fq.gz -o grep_s00_C.bs.bam 2>grep_s00_C.bs.log &
samtools sort -l 9 grep_s00_C.bs.bam grep_s00_C.bs.sort	# UnPaired与病毒比

samtools view -bS grep_s00_C.bs.virsam.gz > grep_s00_C.bs.virsam.bam && \
samtools sort -l 9 grep_s00_C.bs.virsam.bam grep_s00_C.bs.virsam.sort &	# UnPaired与人比

./getfuse.pl 200 grep_s00_C.bs.virsam.sort.bam grep_s00_C.bs.sort.bam
