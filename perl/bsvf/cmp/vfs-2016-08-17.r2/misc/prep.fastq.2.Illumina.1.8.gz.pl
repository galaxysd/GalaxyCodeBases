#!/usr/bin/perl
###########################################################################
##  ViralFusionSeq
##  Software for discovering and annotating viral integration event and
##    fusion transcript
##  
##  Version 1.0 -- November 15, 2012
##  
##  Copyright (C) 2012 by Jing-Woei Li & Raymond Wan, All rights reserved.
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


#This script does 2 things
#(1): Briefly validate Fastq files 
#(2): Convert the fastq to Casava v.1.8 fq and optionally gzip for VFS

##Please use un-compressed fastq files as inputs

#use strict;
use FileHandle;
$| = 1; #disable buffer

my $usage = qq{Usage: $0 <Forward fastq> <Reverse fastq (optional)>
};
die($usage) if (@ARGV < 1);
my $if_gzip = 0;
print STDERR "Do you want to gzip the fastq files?\n";
print STDERR "If your system is short of disk space, select YES. Gzip will save approximate 2/3 of the space required\n";
print STDERR "For maximum VFS performance, select NO\n";
$ans = <STDIN>;
if (($ans =~ /y/gi)||($ans =~ /1/gi)){
	$if_gzip = 1;
	print STDERR "Going to gzip the fastq(s)\n";
} else {
	print STDERR "Will let the fastq(s)un-compressed\n";
}

my $overall_check = 1;
for my $fastq_file (@ARGV) {
	my $status = validate_fastq($fastq_file);
	if ($status != 1){
		$overall_check = 0;
	}
}

my @cmd = ();
if ($overall_check == 1){
	for my $fastq_file (@ARGV) {
		open my $INPUT, '<', $ARGV[0];
		my @input_filename_a = split (/\./, $fastq_file);
		my $range = 1 .. -2;
		my $fastq_file_prefix = join (".", @input_filename_a[$range]);
		my $file_out = $fastq_file_prefix.".vfs_in.fq";
		print "Converting $fastq_file to $file_out\n";
		open my $wFH, '>', $file_out;
		my $i = 1;
		while(<$INPUT>){
			chomp;
			print $wFH '@'.$i." ".'vfs'."\n";
			my $seq = <$INPUT>; print $wFH $seq;
			my $aux = <$INPUT>; print $wFH $aux;
			my $qual = <$INPUT>; print $wFH $qual;
			$i++;
		}
		close ($INPUT);
		close ($wFH);
		if ($if_gzip == 1){
			my $cmd = "gzip -9 -f $file_out";
			print $cmd."\n";
			system($cmd);
		}
	}
} else {
	print STDERR "Fastq failed validation. Please check your fastq files\n";
}

sub validate_fastq {
	my $inputFastq = $_[0];
	print "Validating fastq file: $inputFastq ";
	my $pass = 1;
	my $stage = 1; #1=id, 2=sequence, 3=quality score
	my $line = 0;
	my $validSequences = 0; my $valid_bases; my $valid_bases_qual;
	my $readname;
	my $readlength_max = 1;
	my $readlength_min = 99999;
	my $sequence = "";
	my $qscore = "";
	my $thisqscore;
	my $this_char;
	open my $log_fh, '>', "$inputFastq.validate.log";
	if ((-e $inputFastq)){
		print $log_fh "$inputFastq exists\n";
	}
	else {
		my $error = "$inputFastq does not exists!\n";
		process_error($error, $log_fh);
	}
	if ((-r $inputFastq)){
		print $log_fh "Reading $inputFastq OK\n";
	}
	else {
		my $error = "$inputFastq is not readable!\n";
		process_error($error, $log_fh);
	}
	open my $FH, '<', $inputFastq;
	while(<$FH>) {
		$line++;
		chomp($_);
		if($stage == 1) { # id line
			unless($_ =~ m/^@.+$/){
				print $log_fh "Fail: Invalid format on ID line $line\n";
				$pass = 0;
			}
			$stage = 2;
			undef $readname;
			$readname = substr($_, 1);
			undef $sequence;
			$sequence = "";
			undef $qscore;
			$qscore = "";
		} 
		elsif($stage == 2) { # sequence or + line
			if($_ =~ m/^\+$/) { # + line with optional description
				if($sequence eq ""){ # +, previous line should be sequence
					print $log_fh "Fail: Valid sequence missing before line $line\n"; 
					$pass = 0;
				}
				$stage = 3;
			} 
			elsif($_ =~ m/^\+.+$/) {
				if (/^\+$readname/){
					$stage = 3; #readname after + (optional)
				}
				else {
					print $log_fh "description for $readname: $_"; 
					$stage = 3;
				}
			} 
			elsif($_ =~ m/^[ACTGNURYSWKMBDHV]+$/) { #sequence
				$sequence = $sequence . $_;
				my $base_len = length ($sequence);
				if ($base_len > $readlength_max){
					$readlength_max = $base_len;
				}
				if ($base_len < $readlength_min) {
					$readlength_min = $base_len;
				}
				$valid_bases += $base_len;
			} 
			else { 
				 print $log_fh "Fail: Invalid sequence format on line $line\n";
				 $pass = 0;
			}
		} 
		
		elsif($stage == 3) { # quality score line
			$qscore = $qscore . $_;
			if(length($sequence) == length($qscore)){
				$stage = 1;
				$validSequences++;
			} elsif(length($sequence) < length($qscore)) {
				print $log_fh "Fail: Quality string longer than sequence string on line $line\n";
				$pass = 0;
			}
		}
	}
	unless(($stage == 1) && ($line > 3) && ($validSequences > 0)){
		print $log_fh "Fail: Last entry incomplete\n";
		$pass = 0;
	}
	close($FH);

	if ($pass == 1){
		print $log_fh "Passed: File $inputFastq is a VALID fastq with $validSequences sequences with $valid_bases bases, $valid_bases_qual base quality call and $line lines.\nMax readlength: $readlength_max\nMin readlength: $readlength_min\n";
		print "OK\n";
	} else {
		print STDERR "Error\n";
	}
	close ($log_fh);
	return $pass;
}