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

die "Usage: $0 <single_sam_bam_file> <soap.coverage depthsingle> <output_prefix> <min_overlap> <min_depth>\n" if @ARGV<2;
my $name=shift;
my $cvg=shift;
my $out=shift;
my $overlap=shift;
my $mindepth=shift;
my (%LENGTH,%Stat,%N5090);
if ($name =~ /\.bam$/) {
	open IN,'-|',"$SAMTOOLSBIN view $SAMFILTER $name" or die "Error opening $name : $!\n";
} elsif ($name =~ /\.sam\.gz$/) {
	open IN,'-|',"$SAMTOOLSBIN view $SAMFILTER -S $name" or die "Error opening $name : $!\n";
} elsif ($name =~ /\.sam$/) {
	open IN,'-|',"$SAMTOOLSBIN view $SAMFILTER -S $name" or die "Error opening $name : $!\n";
} else {
	die "[x]Unsupport file type.";
}

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
my $fh=openfile($cvg);
my ($genome,%DepDatChr);
while(<$fh>) {
	s/^>//;
	/^(\S+)/ or next;
	my $SeqName = $1;
	$/=">";
	$genome=<$fh>;
	$/="\n";
	my $i=0;
	while($genome=~/(\d+)/gc) {
		$DepDatChr{$SeqName}->[$i]+=$1 if ($1!=65536);
		++$i;
	}
	$genome='';
	last;
}
close $fh;

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


my ($lastB,$lastE)=(0,0);
while (<IN>) {
	next if /^@\w\w\t\w\w:/;
	chomp;
	my @read1=split /\t/;
	my $cigar=$read1[5];
	next unless $cigar =~ /M/;
	#next unless ($read1[1] & 3) == 3;  # paired + mapped in a proper pair
	next if $read1[1] >= 256;   # not primary || QC failure || optical or PCR duplicate
	#next unless $read1[5] =~ /^(\d+)M$/;
      # http://davetang.org/muse/2011/01/28/perl-and-sam/
      my $position = '0';
      while ($cigar !~ /^$/){
         if ($cigar =~ /^([0-9]+[MIDSH])/){
            my $cigar_part = $1;
            if ($cigar_part =~ /(\d+)M/){
               $position += $1;
            } elsif ($cigar_part =~ /(\d+)[IH]/){
               #$position += $1;
            } elsif ($cigar_part =~ /(\d+)D/){
               $position += $1;
            } elsif ($cigar_part =~ /(\d+)S/){
               die "[!]Not ready for this!\n";
               #my $insertion = 'x' x $1;
               #substr($new_ref,$position,0,$insertion);
               #$position += $1;
            }
            $cigar =~ s/$cigar_part//;
         } else {
            die "Unexpected cigar: $cigar\n";
         }
      }
	my $end=$read1[3]+$position;

	my ($GapFlag,@Gaps)=(0);
	for my $p ($read1[3] .. $end) {
		if ($DepDatChr{$read1[2]}->[$p] < $mindepth) {
			push @Gaps,[$p,-1] unless $GapFlag;
			$GapFlag=1;
		}
		if ($GapFlag and $DepDatChr{$read1[2]}->[$p] >= $mindepth) {
			$Gaps[-1]->[1]=$p;
			$GapFlag=0;
		}
	}
	if ($GapFlag) {
		$Gaps[-1]->[1]=$end+1;
	}
	if ($lastB==0) {	# just the 1st line.
		($lastB,$lastE)=($read1[3],$end);
		if (@Gaps>0) {
			$lastB=$Gaps[0]->[1];
		}
		next;
	}
	if (($lastE-$read1[3] +1 > $overlap) and ($end > $lastE)) {
		;
	}
}
close IN;

