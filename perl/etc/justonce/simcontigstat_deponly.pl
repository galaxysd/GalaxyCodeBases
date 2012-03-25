#!/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 1.0.0 @ 20120309
=cut
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );

my $SAMTOOLSBIN="samtools";
#$SAMTOOLSBIN="/ifs1/ST_ASMB/USER/yuanjy/huxuesong/tmp/group/rev/test/samtools";
my $SAMFILTER="-f 3 -F 1536";
$SAMFILTER="-F 1536";

my $MINCONTIFLEN=100;

die "Usage: $0 <soap.coverage depthsingle> <output_prefix> <min_depth>\n" if @ARGV<2;
my $cvg=shift;
my $out=shift;
my $mindepth=shift;
my (%LENGTH,%Stat,%N5090);

open O,'>',"$out.$mindepth.depthcontig" or die "Error: $!\n";
print O "[Contig]\n#min_depth=$mindepth\n#Min_Contig_Length=$MINCONTIFLEN\n#ChrID\tStart\tEnd\tLen\n";

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

sub don5090($$$) {
	my ($chr,$Len,$href)=@_;
	my $n50 = -1;
	my $n90 = -1;
	my $s = 0;
	my $count=0;
	my @Length = sort {$b <=> $a} keys %{$href};
	for my $l (@Length){
		$s += $l*$$href{$l};
		$count += $$href{$l};
		if($n50 == -1 and $s >= $Len * 0.5){ $n50 = $l; }
		if($n90 == -1 and $s >= $Len * 0.9){ $n90 = $l; last; }
	}
	$N5090{$chr}=[$n50,$n90,$count,$Length[0],$Length[-1]];
	return " [$count,$Length[0],$Length[-1]] $n50,$n90 ";
}
sub dostat($$) {
	my ($chr,$aref)=@_;
	my ($lastB,$lastE)=(0,0);
	my ($GapFlag,@Contig)=(0);
print "\n>$chr\n";
	for my $p (0 .. $#$aref) {
print " $p:$aref->[$p]";
		if ($aref->[$p] >= $mindepth) {
			push @Contig,$p unless $GapFlag;
print "\nB[$p:$aref->[$p]]";
			$GapFlag=1;
		}
		if ($GapFlag and $aref->[$p] < $mindepth) {
			push @Contig,$p;
print "E\n";
			$GapFlag=0;
		}
	}
	push @Contig,$#$aref if $GapFlag;
	my $count=(scalar @Contig)/2;
	while (@Contig) {
		my $a=shift @Contig;
		my $b=shift @Contig;
		my $contiglen = $b-$a+1;
		next if $contiglen < $MINCONTIFLEN;
		++$Stat{'_All_'}{$contiglen};
		++$Stat{$chr}{$contiglen};
		$LENGTH{'_All_'} += $contiglen;
		$LENGTH{$chr} += $contiglen;
		print O join("\t",$chr,$a+1,$b+1,$contiglen),"\n";
	}
	my $n59ret=&don5090($chr,$LENGTH{$chr},$Stat{$chr});
	return $count.$n59ret.$LENGTH{$chr};
}

my $fh=openfile($cvg);
my ($genome,@DepDatChr);
while(<$fh>) {
	s/^>//;
	/^(\S+)/ or next;
	my $SeqName = $1;
	@DepDatChr=();
	print STDERR "$SeqName: ";
	$/=">";
	$genome=<$fh>;
	$/="\n";
	my $i=0;
	while($genome=~/(\d+)/gc) {
		$DepDatChr[$i]+=$1 if ($1!=65536);
		++$i;
	}
	my $cnt=&dostat($SeqName,\@DepDatChr);
	$genome='';
	print STDERR "$cnt\n";
}
close $fh;

print O "\n\n[Stat]\n#ChrID\tLen\tCount\n";
for my $chr (sort keys %Stat) {
	for my $len (sort {$b<=>$a} keys %{$Stat{$chr}}) {
		print O join("\t",$chr,$len,$Stat{$chr}{$len}),"\n";
	}
}

my $chr = '_All_';
don5090($chr,$LENGTH{$chr},$Stat{$chr});
print O "\n\n[N50_N90]\n#ChrID\tN50\tN90\tCount\tMaxLen\tMinLen\tTotal_Length\n";
for my $chr (sort keys %N5090) {
	print O join("\t",$chr,@{$N5090{$chr}},$LENGTH{$chr}),"\n";
}
close O;
