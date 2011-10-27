#!/usr/bin/perl
use strict;
use Getopt::Long;
my %opts;

my $usage=<<USAGE; 

Function: (1) split scaffold sequences into contig sequnences;
	(2) get contig coordinates and length;  
	(3) get gap coordinates and length

Contact: Fan wei,<fanw\@genomics.org.cn>

Version: 1.0, release 2010-1-7

Usage: perl split_scaffold.pl [options] <input_scaffold_file.fa>
	-N <num>    number of consecutive Ns for split site,default N=1
	-out <str>  three alternative types: contig_seq, contig_coor, gap_coor, default=contig_seq
	-cutoff <num>  set the size cutoff for contig, only contig > cutoff will be output
	-detail		output middle information to screen
	-help		output help information to screen

Example: perl split_scaffold.pl scaffold.fa > contig.fa
Note that the -cutoff option can only be co-used with contig_seq, contig_coor.

USAGE

GetOptions(\%opts,"N:i", "out:s", "cutoff:i", "detail!","help!");
die $usage if ( @ARGV==0 || defined($opts{"help"}));

####################################################
################# Main  Function ###################
####################################################
my $file=shift;
my $cut_N= (exists $opts{N}) ? $opts{N} : 1;
my $cutoff =(exists $opts{cutoff}) ? $opts{cutoff} : 0;
my $out_type = (exists $opts{out}) ? $opts{out} : "contig_seq";

open(IN, $file) || die ("can not open $file\n");
$/=">"; <IN>; $/="\n";
while (<IN>) {
	my $chr=$1 if(/^(\S+)/);
	$/=">";
	my $seq=<IN>;
	chomp $seq;
	$seq=~s/\s//g;
	$seq = uc $seq;
	$/="\n";
	
	print STDERR "\nSplit $chr\n" if (exists $opts{detail});
	
	my @frag;
	my ($start,$end,$frag,$pos);
	
	#Remove Ns at heads
	while ($seq=~s/^(N+)//g) {
		$pos+=length($1);
	}
	
	##cut scaffold to contigs
	while ($seq=~s/^([^N]+)(N*)//) {
		my $len_A=length($1);
		my $len_N=length($2);
		$pos += $len_A+$len_N;
		if ($len_N >= $cut_N || !$seq) {
			$frag .= $1;
			$end=$pos-$len_N;
			$start = $end - $len_A + 1;
			my $gap_end = $pos;
			my $gap_start = $pos - $len_N + 1;
			push @frag,[$frag, $start,$end, $len_A, $gap_start, $gap_end, $len_N];
			$frag="";
		}else{
			$frag .= $1.$2;
		}
	
	}
	
	##output the result in various format
	my $frag_num = @frag;
	for (my $i=0; $i<$frag_num; $i++) {
		Display_seq(\$frag[$i][0]);
		my $mark = $i+1;
		my $contig_id = "$chr\_$mark";
		
		print ">$contig_id  $chr $frag[$i][1] $frag[$i][2] $frag[$i][3]\n$frag[$i][0]" if ($out_type eq "contig_seq" && $frag[$i][3]>=$cutoff);
		print "$chr $frag[$i][1] $frag[$i][2] $frag[$i][3]\n" if ($out_type eq "contig_coor"  && $frag[$i][3]>=$cutoff);
		print "$chr $frag[$i][4] $frag[$i][5] $frag[$i][6]\n" if ($out_type eq "gap_coor"  && $frag[$i][6]>0);
	}

	
}
close(IN);




#display a sequence in specified number on each line
#usage: disp_seq(\$string,$num_line);
#		disp_seq(\$string);
#############################################
sub Display_seq{
	my $seq_p=shift;
	my $num_line=(@_) ? shift : 50; ##set the number of charcters in each line
	my $disp;

	$$seq_p =~ s/\s//g;
	for (my $i=0; $i<length($$seq_p); $i+=$num_line) {
		$disp .= substr($$seq_p,$i,$num_line)."\n";
	}
	$$seq_p = ($disp) ?  $disp : "\n";
}
#############################################
