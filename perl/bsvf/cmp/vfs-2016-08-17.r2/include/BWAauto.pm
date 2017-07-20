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

package BWAauto;

use diagnostics;
use strict;
use warnings;

$| = 1; #disable buffer

use FileHandle;
use Exporter;
use Statistics::Descriptive;

use FindBin;
use lib "$FindBin::Bin/include";
use General; # Module for general subroutines

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

our $VERSION     = 1.00;
our @ISA         = qw(Exporter);
our @EXPORT      = qw(bwa_pipeline);
our @EXPORT_OK   = qw(bwa_pipeline);
our %EXPORT_TAGS = ( DEFAULT => [qw(&bwa_pipeline)]);

#### default input ####
my $remove_temp = 0;
my $illumina13Q = 0;
my $bwa_illumina13 = ' ';
if ($illumina13Q == 1){
	$bwa_illumina13 = ' -I ';
}
my $backtrim = 0; #Assume *.fastq as suffix of input fastq
my $obtain_IS = 0; #Do not return insert size by default
my $bwa_trimQ = 0;
my $fast_run = 1;
#########################
my @log_files_generated = ();

sub bwa_pipeline {
	my $curr_func = (caller(0))[3];
	my $ref_genome = $_[0];
	my $SAM_prefix = $_[1];
	my $threads = $_[2];
	my $F_fq = $_[3];
	my $R_fq = $_[4];
	my $bwa = $_[5];
	my $obtain_IS = $_[6] if (defined $_[6]);
	my $samtools = $_[7];
	my $infer_Isize_median = "NA";
	my $Rinfer_Isize = "NA";
	if (defined $R_fq){ # PE mode
		#aln
		print STDERR "Running BWA in paired-ended mode\n";
		print STDERR "Aligning $F_fq\n";
		my $outputSAI_F_tmp = backtrim_string($F_fq, $backtrim);
		my $outputSAI_F = "$outputSAI_F_tmp.sai";
		my $F_aln_log = "$SAM_prefix.Faln.stderr";
		push (@log_files_generated, $F_aln_log);
		my @Sampe_Faln_cmd = ($bwa, "aln", "-q", "$bwa_trimQ"."$bwa_illumina13", "-t", $threads, $ref_genome, $F_fq, ">$outputSAI_F", "2>$F_aln_log");
		print join (" ", @Sampe_Faln_cmd), "\n";
		ExecCmd (\@Sampe_Faln_cmd, $DEBUG_MODE);

		print STDERR "Aligning $R_fq\n";
		my $outputSAI_R_tmp = backtrim_string($R_fq, $backtrim);
		my $outputSAI_R = "$outputSAI_R_tmp.sai";
		my $R_aln_log = "$SAM_prefix.Raln.stderr";
		push (@log_files_generated, $R_aln_log);
		my @Sampe_Raln_cmd = ($bwa, "aln", "-q", "$bwa_trimQ"."$bwa_illumina13", "-t", $threads, $ref_genome, $R_fq, ">$outputSAI_R", "2>$R_aln_log");
		print join (" ", @Sampe_Raln_cmd), "\n";
		ExecCmd (\@Sampe_Raln_cmd, $DEBUG_MODE);
		#aln

		#sampe
		print STDERR "BWA sampe\n";
		my $sampe_log = "$SAM_prefix.sampe.stderr";
		push (@log_files_generated, $sampe_log);
		my @sampe_cmd;
		if ($obtain_IS == 1){
			if ($fast_run == 1){
				@sampe_cmd = split (/ /, "$bwa sampe -P $ref_genome $outputSAI_F $outputSAI_R $F_fq $R_fq >/dev/null 2>$sampe_log");
			} else {
				@sampe_cmd = split (/ /, "$bwa sampe -P $ref_genome $outputSAI_F $outputSAI_R $F_fq $R_fq >$SAM_prefix.sam 2>$sampe_log");
			}
		} else {
			@sampe_cmd = split (/ /, "$bwa sampe -P $ref_genome $outputSAI_F $outputSAI_R $F_fq $R_fq | $samtools view -Sbu - | $samtools sort -m 2000000000 - $SAM_prefix.sorted 2>$sampe_log");
		}
		print join (" ", @sampe_cmd), "\n";
		ExecCmd (\@sampe_cmd, $DEBUG_MODE);

		print STDERR "Cleaning up Mapping files\n";
		#cleanup SAI
		unlink ($outputSAI_F);
		unlink ($outputSAI_R);
		#print "OK\n";
		#obtain insert size
		if ($obtain_IS == 1){
			$Rinfer_Isize = get_median_ext_Isize($SAM_prefix, "$SAM_prefix.sampe.stderr");
		}
	}
	else { #single-ended mode
		print STDERR "Running BWA in single-ended mode\n";
		print STDERR "Aligning $F_fq\n";
		my $outputSAI_tmp = backtrim_string($F_fq, $backtrim);
		my $outputSAI_F = "$outputSAI_tmp.sai";
		my $log = "$SAM_prefix.Faln.stderr";
		push (@log_files_generated, $log);
		my @cmd_samse_aln = ($bwa, "aln", "-t", $threads, "-q", "$bwa_trimQ"."$bwa_illumina13"."$ref_genome", $F_fq, ">$outputSAI_F", "2>$log");
		print join (" ", @cmd_samse_aln), "\n";
		ExecCmd (\@cmd_samse_aln, $DEBUG_MODE);

		#samse
		my $samse_log = "$SAM_prefix.samse.stderr";
		push (@log_files_generated, $samse_log);
		my @samse_cmd = "bwa samse $ref_genome $outputSAI_F $F_fq | $samtools view -Sb - | $samtools sort -m 2000000000 - $SAM_prefix.sorted 2>$samse_log";
		print join (" ", @samse_cmd), "\n";
		ExecCmd (\@samse_cmd, $DEBUG_MODE);

		#cleanup SAI
		print STDERR "Cleaning up Mapping files\n";
		unlink ($outputSAI_F);
	}

	#invoke samtools
	if ($obtain_IS == 1){
		unless ($fast_run == 1){
			my @cmd;

			unlink ("$SAM_prefix.sam");

			@cmd = split (/ /, "$samtools view -Sb $SAM_prefix.sam >$SAM_prefix.bam");
			ExecCmd (\@cmd, $DEBUG_MODE);
			unlink ("$SAM_prefix.sam");

			@cmd = split (/ /, "$samtools sort -m 2000000000 $SAM_prefix.bam $SAM_prefix.sorted");
			ExecCmd (\@cmd, $DEBUG_MODE);
			unlink ("$SAM_prefix.bam");

			@cmd = split (/ /, "$samtools flagstat $SAM_prefix.sorted.bam >$SAM_prefix.sorted.bam.flagstat");
			ExecCmd (\@cmd, $DEBUG_MODE);

			@cmd = split (/ /, "$samtools index $SAM_prefix.sorted.bam");
			ExecCmd (\@cmd, $DEBUG_MODE);
		}
	}
	return $Rinfer_Isize if ($obtain_IS == 1);
}

#clean up log files upon request
if ($remove_temp == 1){
	for (@log_files_generated){
		unlink ($_);
	}
}

# Function
sub get_median_ext_Isize {
	my $SAM_prefix = $_[0];
	my $file = $_[1];
	my @Isizes = ();
	my @Isizes_stdev = ();
	open my $Isizefh, '>', "$SAM_prefix.isize" or die "Can't open $SAM_prefix.isize for writing : $!";
	open my $fh, '<', $file or die "Can't open $file for reading : $!";
	while (<$fh>){
		chomp;
		if (/infer_isize/){
			if (/inferred external isize/){
				$_ =~ m/(\d+\.\d+)\s\+\/-\s(\d+\.\d+)/gi;
				my ($Isize, $stdev) = ($1, $2);
				print $Isizefh "$Isize\t$stdev\n";
				push (@Isizes, $Isize);
			}
		}
	}
	close ($fh);

	#obtain median of inferred insert size
	my $fullstat = Statistics::Descriptive::Full->new();
	$fullstat->add_data(@Isizes);
	my $median = $fullstat->median();
	my $Rinfer_Isize = sprintf("%.0f", $median);
	print $Isizefh "inferred insert size: $median\n";
	print $Isizefh "rounded insert size: $Rinfer_Isize\n";
	close ($Isizefh);
	return $Rinfer_Isize;
}

sub round{
    my($number) = shift;
    return int($number + .5 * ($number <=> 0));
}

sub backtrim_string{
	my @input = @_;
	my $inputstring = $input[0];
	my $offset = $input[1];
	my $length = length ($inputstring);
	my $span = $length - $offset;
	my $out = substr $inputstring, 0, $span;
	return $out;
}

1;
