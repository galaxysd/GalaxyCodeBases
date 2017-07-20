#!/usr/bin/perl
###########################################################################
##  ViralFusionSeq
##  Software for discovering and annotating viral integration event and
##    fusion transcript
##  
##  Version 2.0
##  
##  Copyright (C) 2016 by Jing-Woei Li & Raymond Wan, All rights reserved.
##  Contact:  marcoli@cuhk.edu.hk, rwan@cuhk.edu.hk
##  Organization:  Hong Kong Bioinformatics Centre, School of Life Sciences, The
##                 Chinese University of Hong Kong, Shatin, NT,
##                 Hong Kong SAR
##  
##  This file is part of ViralFusionSeq.
##  
##  ViralFusionSeq is free software; you can redistribute it and/or 
##  modify it under the terms of the GNU General Public License 
##  as published by the Free Software Foundation; either version 
##  3 of the License, or (at your option) any later version.
##  
##  ViralFusionSeq is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##  
##  You should have received a copy of the GNU General Public 
##  License along with ViralFusionSeq; if not, see 
##  <http://www.gnu.org/licenses/>.
###########################################################################

##  $LastChangedDate: 16th Aug 2016 $
##  $LastChangedRevision: 1289 $

package ReadPreprocess;

#Dynamic reads trimming
#for illumina pipeline v.1.8 reads

##  Modules included with this program
use General; # Module for general subroutines
use diagnostics;
use strict;
use warnings;

use Exporter;
use FileHandle;
use PerlIO::gzip;

$| = 1; #disable buffer
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = qw(preprocess_fastq run_dynamic_trim_fastq split_n_take_first union_filter dynamic_trim_fastq_tag process_error usage);
@EXPORT_OK   = qw(preprocess_fastq run_dynamic_trim_fastq split_n_take_first union_filter dynamic_trim_fastq_tag process_error usage);
%EXPORT_TAGS = ( DEFAULT => [qw(&preprocess_fastq &run_dynamic_trim_fastq &split_n_take_first &union_filter &dynamic_trim_fastq_tag &process_error &usage)]);

sub preprocess_fastq { #(#file_in, expected sequencing length, phred scale for reduction (Sanger 33, Solexa 64), read length emit threshold)
	my $curr_process =(caller(0))[3];
	print STDERR "executing $curr_process.\n";
	my $file_in_one = $_[0];
	my $file_in_two = $_[1];
	my $phred_scale = $_[2];
	my $desired_quality = $_[3];
	my $emit_threshold = $_[4];  #minimum read length after trimming
	my $PE = 1;
	if(!defined($file_in_two)){
		$PE = 0;
		print STDERR "Reverse read is not defined. VFS will directly process $file_in_one\n";
			print STDERR "Fastq file: $file_in_one"."\n"."Phred-Scale: $phred_scale"."\n"."Desired Base call Q: $desired_quality"."\n";
	} else {
			print STDERR "Fastq file: $file_in_one and $file_in_two"."\n"."Phred-Scale: $phred_scale"."\n"."Desired Base call Q: $desired_quality"."\n";
	}

	if ((!defined($phred_scale)) || (($phred_scale != 64) && ($phred_scale != 33))){
		if (!defined($phred_scale)){
			print STDERR "Phred Scale is NOT defined\n";
		}
		if (($phred_scale != 64) && ($phred_scale != 33)){
			print STDERR "Phred scale must be Numerical value and has to be(64/33)\n";
		}
		usage();
		exit;
	}
	if((!defined($desired_quality)) || ($desired_quality !~ m/^[0-9]+$/)){
		print STDERR "Desired base call quality must be integer";
		usage();
		exit;
	}
	my @out_files = ();
	if ($PE == 1){
		@out_files = @{run_dynamic_trim_fastq_PE($file_in_one, $file_in_two, $phred_scale, $desired_quality, $emit_threshold)};
	} else {
		@out_files = @{run_dynamic_trim_fastq_SE($file_in_one, $phred_scale, $desired_quality, $emit_threshold)};
	}
	return \@out_files;
}

sub run_dynamic_trim_fastq_SE {
	my $file_in_one = $_[0];
	my $phred_scale = $_[1];
	my $desired_quality = $_[2];
	my $emit_threshold = $_[3];
	my @out_files = ();
	dynamic_trim_fastq_tag($file_in_one, $phred_scale, $desired_quality, $emit_threshold);
	my $out_prefix = split_n_take_first($file_in_one);
	my $F_output_file = "$out_prefix.dynamic.min.$emit_threshold.fq.gz";
	rename "$file_in_one.trim.temp.gz", $F_output_file;
	push(@out_files, $F_output_file);
	return \@out_files;
}

sub run_dynamic_trim_fastq_PE {
	my $file_in_one = $_[0];
	my $file_in_two = $_[1];
	my $phred_scale = $_[2];
	my $desired_quality = $_[3];
	my $emit_threshold = $_[4];
	my @out_files = ();
	my @one_discard_name = @{dynamic_trim_fastq_tag($file_in_one, $phred_scale, $desired_quality, $emit_threshold)};
	my @two_discard_name = @{dynamic_trim_fastq_tag($file_in_two, $phred_scale, $desired_quality, $emit_threshold)};
	my %union =();
	foreach my $ele (@one_discard_name) {
		$union{$ele} = 1;
	}
	foreach my $ele (@two_discard_name) {
		$union{$ele} = 1;
	}
	#revisit each file to synchronize the reads on two ends
	my $F_output_file = union_filter("$file_in_one.trim.temp.gz", \%union, 'dynamic',1, $emit_threshold);
	unlink ("$file_in_one.trim.temp.gz");
	my $R_output_file = union_filter("$file_in_two.trim.temp.gz", \%union, 'dynamic',1, $emit_threshold);
	unlink ("$file_in_two.trim.temp.gz");
	push(@out_files, $F_output_file, $R_output_file);
	return \@out_files;
}

sub split_n_take_first {
	my $in = $_[0];
	my @splitted = split (/\./, $in);
	my $rtn = $splitted[0];
	return $rtn;
}

##  Input and output files are both gzipped
sub union_filter {
	my $file_in = $_[0];
	print STDERR "Union Filtering $file_in\t";
	my %union = %{$_[1]};
	my $description = $_[2];
	my $trim = $_[3];
	my $trim_emit_threshold = $_[4];
	my $out_prefix = "";
	if ($trim){
		$out_prefix = split_n_take_first($file_in);
	}
	my $readname;  my $seq;  my $plus;  my $quality_score; #static info

	my $output_file = "$out_prefix.$description.min$trim_emit_threshold.fq.gz";

	open (my $gz_out, ">:gzip", "$output_file") or die "Can't open $output_file (gzipped) for writing : $!";
	open (my $gz_in, "<:gzip", "$file_in") or die "Can't open $file_in (gzipped) for reading : $!";

	while (<$gz_in>) {
		$readname = $_;
		chomp ($readname);
		$seq = <$gz_in>;  # read next line, which is the sequence
		$plus = <$gz_in>; # plus sign plus any description
		$quality_score = <$gz_in>;  # read quality scores
		#my $read = substr($readname, 0, -2); for xxxx/1, xxxx/2
		my @sr = split (/ /, $readname); #@HWI-ST993:188:C0CA8ACXX:6:1101:1847:1913 1:N:0:ACTTGA "@HWI-ST993:188:C0CA8ACXX:6:1101:1847:1913" is extracted
		my $read = $sr[0];
		unless (defined $union{$read}){
			my $outline = $readname."\n".$seq.$plus.$quality_score;
			print $gz_out $outline;
		}
	}

	close ($gz_in);
	close ($gz_out);
	print STDERR "OK\n";

	#return prefix name
	return $output_file;
}

sub dynamic_trim_fastq_tag { #(File_in, phred_scale, desiredQ, emit_threshold)
	print STDERR "Dynamic tagging reads ";
	my $file_in = $_[0];
	my $phred_scale = $_[1];
	my $desQ = $_[2];
	my $emit_threshold = $_[3];
	my $readname;  my $seq;  my $plus;  my $quality_score; my $readlength; #static info
	my $pos;  my $maxPos;  my $area;  my $maxArea; my $substr_len;
	my $total_processed_reads; my $discard_read_number;
	my @discard_read_names = ();

	open (my $gz_out, ">:gzip", "$file_in.trim.temp.gz") or die "Can't open $file_in.trim.temp.gz (gzipped) for writing : $!";

	my $FH;
	my $is_gzip_input = DetermineFileType ($file_in);

	if ($is_gzip_input) {
		open ($FH, "<:gzip", $file_in) or die "Can't open $file_in (gzipped) for reading : $!";
	}
	else {
		open $FH, '<', $file_in or die "Can't open $file_in for reading : $!";
	}

	while (1) {
		$readname = <$FH>;
		if (!defined $readname) {
			last;  ##  End of file
		}
		$seq = <$FH>;
		$plus = <$FH>;
		$quality_score = <$FH>;

		chomp $readname;
		chomp $seq;
		chomp $plus;
		chomp $quality_score;

		if (length($seq) == length ($quality_score)){
			$readlength = length($quality_score);
		}
		else {
			die('Sequence length not equal to quality length. Validate the FastQ file first');
		}
		$pos = $readlength;
		$maxPos = $pos;
		$area = 0;
		$maxArea = 0;
		while ($pos>0 and $area>=0) {
			my $no = (ord(substr($quality_score,$pos-1,1))- $phred_scale);
			$area += $desQ - (ord(substr($quality_score,$pos-1,1))- $phred_scale);
			if ($area > $maxArea) {
				$maxArea = $area;
				$maxPos = $pos;
			}
			$pos--;
		}
		if ($maxPos == $readlength){
			$substr_len = $maxPos;
		}
		else {
			$substr_len = $maxPos-1;
		}
		#trim from start of read to position before position where area reached a maximum
		my $ori_seq = $seq;
		my $ori_quality_score = $quality_score;
		$seq = substr($seq,0,$substr_len);
		$quality_score = substr($quality_score,0,$substr_len);

		unless (length($seq) < $emit_threshold){ #have checked previously read length == quality length
			my $out_line = $readname."\n".$seq."\n".$plus."\n".$quality_score."\n";
			print $gz_out $out_line;
		}
		else {
			my @sr = split (/ /, $readname); #this part should be modified in case sequence reads were generated from illumina pipeline v.1.3. ReadID naming issue.
			push (@discard_read_names, $sr[0]);
		}
	}

	close ($gz_out);
	close ($FH);

	print STDERR "OK\n";
	return \@discard_read_names;
}

sub process_error {
	my $msg = $_[0];
	my $log_fh = $_[1];
	print $log_fh $msg;
	warn $msg;
	exit;
}

sub usage {
	print qq{Usage: $0 <F_readfile>, <R_readfile (optional)>, <phredQ>, <desiredQ>, <emitThreshold>
	};
}

1;
