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

#ddx \%ChrLen; ddx \%Genome;
ddx $finHum;
ddx $finHBV;
my ($MaxReadLen) = sort {$b<=>$a} ($$finHum[3],$$finHum[4],$$finHBV[3],$$finHBV[4]);

warn "[!]Ref: h=[$$finHum[0]],v=[$$finHBV[0]]. MaxReadLen=[$MaxReadLen]\n";
my $SearchPanRange = 2*$isize + $MaxReadLen;	# incase reads with soft clipping.
=pod
也需要统计reads断点的累计来推测断点。
=cut
sub doPileUp() {
	my ($t)=@_;
	my (%PileUp,%Depth);
}
my $doPanDAT = ['',-$SearchPanRange];	# Chr, RightPos
sub doPan($$$$) {
	my ($ref,$pos,$mappedLen,$id)=@_;
	#my (%PileUp,%Depth);
	if (($ref ne $$doPanDAT[0]) or ($pos > 1 + $$doPanDAT[1] + $SearchPanRange)) {	# Next Pan is coming, analyse last Pan now.
		doPileUp() if $$doPanDAT[0] ne '';
		$$doPanDAT[0] = $ref; $$doPanDAT[1] = $pos+$mappedLen-1;
	}
}
=pod
open INHUM,"-|","samtools view $inhum" or die "Error opening $inhum: $!\n";
while (<INHUM>) {
	chomp;
	my ($id, $flag, $ref, $pos, $mapq, $CIGAR, $mref, $mpos, $isize, $seq, $qual, @OPT) = split /\t/;
	my $mappedSeq = parseCIGAR($seq,$CIGAR);
	my $mappedQual = parseCIGAR($qual,$CIGAR);
	my $mappedLen = length $mappedSeq;
	doPan($ref,$pos,$mappedLen,$id,$mappedSeq,$mappedQual);
}
close INHUM;
=cut
open INVIR,"-|","samtools view -F 4 $invir" or die "Error(1) opening $invir: $!\n";
my (%PileUpVIR,%DepthVIR);
while (<INVIR>) {
	chomp;
	my ($id, $flag, $ref, $pos, $mapq, $CIGAR, $mref, $mpos, $isize, $seq, $qual, @OPT) = split /\t/;
	my $mappedSeq = parseCIGAR($seq,$CIGAR);
	my $mappedQual = parseCIGAR($qual,$CIGAR);
	my $mappedLen = length $mappedSeq;
	if ( $DEBUG > 1 ) {
		doPan($ref,$pos,$mappedLen,$id);
		push @{$PileUpVIR{$ref}},[$pos,$mappedLen,$id,$mappedSeq,$mappedQual,$seq,$CIGAR];
warn "$CIGAR,$ref\t$id,$seq\n$mappedSeq\n";
		my $theref = substr($Genome{$ref},$pos-1,$mappedLen);
		# 用bisofite处理后，含甲基化的C不会变成T，不含甲基化的C变成T。由于他只能吧C变成T，所以，正链就是C->T，而他的反义互补连则是负连的C->T，反映到正链上就是G->A。
		# 比如CG，如果他在watson连上并且没有甲基化，那么bisofte处理后就变成TG，如果是在crick链上变成CA。
		my $therefMet = $theref; $therefMet =~ s/C/T/ig;
		#my $therefRCMet = revcom($theref); $therefRCMet =~ s/C/T/ig; $therefRCMet = revcom($therefRCMet);
		my $therefRCMet = $theref; $therefRCMet =~ s/G/A/ig;
warn "$therefMet\n$therefRCMet\n$theref\n\n";
	} else {
		push @{$PileUpVIR{$ref}},[$pos,$mappedLen,$id,$mappedSeq,$mappedQual];
	}
}
close INVIR;
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

./getfuse.pl 200 grep_s00_C.bs.virsam.sort.bam grep_s00_C.bs.sort.bam >t.log 2>t.err

------

./bwameth.py --calmd -t 24 -p s00_C.bshum.calmd --read-group s00_C --reference HomoGRCh38.fa s00_C.bs_1.fq.gz s00_C.bs_2.fq.gz 2>s00_C.bshum.calmd.err
./bwameth.py --calmd -t 24 -p s00_C.bshbv.calmd --read-group s00_C --reference HBV.AJ507799.2.fa s00_C.bs_1.fq.gz s00_C.bs_2.fq.gz 2>s00_C.bshbv.calmd.err

------

./bwameth.py -t 24 -p s00_C.bshum --read-group s00_C --reference HomoGRCh38.fa s00_C.bs_1.fq.gz s00_C.bs_2.fq.gz 2>s00_C.bshum.err
samtools sort -m 2415919104 -n s00_C.bshum.bam -O bam -T s00_C.bshum.sid >s00_C.bshum.sid.bam 2>tt.err

../getUnPaired.pl s00_C.bshum.sid.bam n3_grep
samtools view -bS n3_grep.vircandi.sam.gz >n3_grep.vircandi.bam
samtools sort -m 2415919104 n3_grep.vircandi.bam -O bam -T n3_grep.vircandi.sort >n3_grep.vircandi.sort.bam

./bwameth.py -t 24 -p n3_grep.vircandi.bshbv --read-group s00_C_grep --reference HBV.AJ507799.2.fa n3_grep.1.fq.gz n3_grep.2.fq.gz 2>n3_grep.vircandi.bshbv.err
samtools sort -m 2415919104 -n n3_grep.vircandi.bshbv.bam -O bam -T n3_grep.vircandi.bshbv.sid >n3_grep.vircandi.bshbv.sid.bam




./bwameth.py -t 24 -p s00_C.bshbv --read-group s00_C --reference HBV.AJ507799.2.fa s00_C.bs_1.fq.gz s00_C.bs_2.fq.gz 2>s00_C.bshbv.err
samtools sort -m 2415919104 -n s00_C.bshbv.bam -O bam -T s00_C.bshbv.sid >s00_C.bshbv.sid.bam

ls -lh s00_C.bshbv.bam s00_C.bshbv.sid.bam s00_C.bshum.bam s00_C.bshum.sid.bam

