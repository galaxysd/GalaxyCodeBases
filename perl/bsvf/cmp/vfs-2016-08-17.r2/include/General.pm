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

package General;

use diagnostics;
use strict;
use warnings;

$| = 1; #disable buffer

##  Modules not provided with this program
use File::Which;
use PerlIO::gzip;  ##  Ubuntu:  libperlio-gzip-perl

use vars qw ($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
our $VERSION     = 1.00;
our @ISA         = qw (Exporter);
our @EXPORT      = qw ($DEBUG_MODE CreateLogDirectory CleanLogDirectory DisplayTime ExecCmd ExpandProgramPath DetermineFileType OpenFastqFileRead GetNextFastqLine GetNextFastqRecord CloseFastqFile CheckPositiveInt CheckPhredQ CheckPath CheckFile CheckProgram CountORFs GetNameORFs GetStartORFs GetEndORFs CheckORF in_range);
our @EXPORT_OK   = qw ($DEBUG_MODE CreateLogDirectory CleanLogDirectory DisplayTime ExecCmd ExpandProgramPath DetermineFileType OpenFastqFileRead GetNextFastqLine GetNextFastqRecord CloseFastqFile CheckPositiveInt CheckPhredQ CheckPath CheckFile CheckProgram CountORFs GetNameORFs GetStartORFs GetEndORFs CheckORF in_range);
our %EXPORT_TAGS = ( DEFAULT => [qw(&CreateLogDirectory &CleanLogDirectory &DisplayTime &ExecCmd &ExpandProgramPath &DetermineFileType &OpenFastqFileRead &GetNextFastqLine &GetNextFastqRecord &CloseFastqFile &&CheckPositiveInt &CheckPhredQ &CheckPath &CheckFile &CheckProgram &CountORFs &GetNameORFs &GetStartORFs &GetEndORFs &CheckORF &in_range)]);


##  Debugging variable; set to 0 or 1 if you want non-debug mode or debug mode
##    Debug mode prints more information while running
our $DEBUG_MODE = 1;

sub DisplayTime {
	 my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
	 my @weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
	 my ($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
	 my $year = 1900 + $yearOffset;
	 my $time = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
	 print STDERR "$time\n";
}

sub in_range {
	my @boundary = @{$_[0]};
	my $start = $_[1];
	my $end = $_[2];
	my $status = 0;
	@boundary = sort {$a <=> $b} @boundary; #min: $array[0]\nmax: $array[-1]
	if (($start >= $boundary[0]) && ($end <= $boundary[-1])){
		$status = 1;
	}
	#print "start: $start\tend: $end\tMax: $boundary[0]\tMin: $boundary[-1]\t$status\n";
	return $status;
}


# Takes two arguments:  an array and an integral flag.
#   Flag's meaning:  1 means verbose, #2 means mock action
sub ExecCmd {
	my $cmd = $_[0];
	my $switch = $_[1];
	if ($switch == 1){
		print STDERR "Executing:  ";
		print STDERR join (" ", @$cmd);
		print STDERR "\n";
	}
	unless ($switch == 2){
		##  We need to execute in scalar context since we occasionally capture the STDERR or STDOUT.
		##  List context is more efficient, but the shell is not called, making STDERR/STDOUT unavailable.
		system (join (" ", @$cmd));
	}
}



##  Expand a program's path if the path given is neither absolute nor relative.
##    If absolute or relative, then return as-is.  If program cannot be found,
##    then print a warning but return it as-is and let the calling program raise
##    an error!
sub ExpandProgramPath {
	my $path = shift;

	if ($path =~ /^\//) {
		##  Absolute path to program; return unchanged
		return $path;
	}
	elsif ($path =~ /^\.\//) {
		##  Program in current directory; return unchanged
		return $path;
	}
	elsif ($path =~ /^\.\.\//) {
		##  Program in previous directory; return unchanged
		return $path;
	}

	my $new_path = which ($path);
	if (!defined $new_path) {
		printf STDERR "WW\tThe program \'%s\' could not be found in the path.\n", $path;
		$new_path = $path;
	}

	return ($new_path);
}


############################################################
#  General functions for accessing FASTQ (or FASTA) records
#    Abstracts away whether the FASTQ file is
#    compressed or not.  Only functions that operate on
#    records cannot work on FASTA data.
############################################################

##  Returns 1 for a gzip file; 0 otherwise
sub DetermineFileType {
	my $filename = shift;

	if ($filename =~ /\.gz$/) {
		return 1;
	}

	return 0;
}


sub OpenFastqFileRead {
	my ($filename, $is_gzip_file) = @_;
	my $fh;

	if ($is_gzip_file) {
		open ($fh, "<:gzip", $filename) or die "Can't open $filename (gzipped) for reading : $!";
	}
	else {
		open ($fh, '<', $filename) or die "Can't open $filename for reading : $!";
	}

	return ($fh);
}


##  Get the next line from a file; abstracts away whether or not it is compressed.
sub GetNextFastqLine {
	my ($fh) = @_;
	my $line;

	$line = <$fh>;
	if (!defined ($line)) {
		return "";
	}

	return ($line);
}


##  Get the next FASTQ record, with fields new line-separated.
##    Returns an empty string if EOF reached.
sub GetNextFastqRecord {
	my ($fh) = @_;
	my $id;
	my $seq;
	my $comm;
	my $qual;

	$id = <$fh>;
	if (!defined ($id)) {
		return "";
	}
	$seq = <$fh>;
	$comm = <$fh>;
	$qual = <$fh>;

	my $record = $id.$seq.$comm.$qual;

	return ($record);
}

sub CloseFastqFile {
	my ($fh) = @_;

	close ($fh);

	return;
}


##  Function is unused; removed by rwan on 2013/02/27
# sub WriteFastqRecord {
# 	my ($fh, $id, $seq, $comm, $qual) = @_;
# 	my $outline = "";
#
# 	$outline = $id."\n".$seq."\n".$comm."\n".$qual."\n";
# 	print $fh $outline;
#
# 	return;
# }


###################################################
#  General functions for logging
###################################################

sub CreateLogDirectory {
	my ($logdir, $debug) = @_;

	if (-e $logdir) {
		print STDERR "Log directory \"$logdir\" already present, cleaning it.\n";
		CleanLogDirectory ($logdir, $debug);
	}
	mkdir $logdir;
	return;
}

sub CleanLogDirectory {
	my ($logdir, $debug) = @_;
	my $files_deleted = unlink glob "$logdir/*";
	if ($debug) {
		#printf "Number of log files deleted:  %u\n", $files_deleted;
	}
	return;
}

###################################################
#  General functions for error handling
###################################################

sub CheckPositiveInt {
	my ($var, $value) = @_;

	if ($value <= 0) {
		return 0;
	}
	return 1;
}

##  Checks phred quality score. 33/64, otherwise use assume NA
sub CheckPhredQ {
	my ($var, $value) = @_;
	my $t = lc($value);
	if ($t =~ /[a-z]/){
		return 1;
	}
	else {
		my $int = int($value);
		if (($int != 33) && ($int != 64)) {
			return 0;
		}
	}
}


##  Check that a path (directory) has been given
sub CheckPath {
	my ($path, $srcname, $subline, $debug) = @_;

	if (!-e $path) {
		if ($debug) {
			printf STDERR "Debug:  Called from \'%s\' at line \'%u\'...\n", $srcname, $subline;
		}
		printf STDERR "Error:  The path \'%s\' does not exist!  Stopped\n", $path;
		return 0;
	}

	if (!-d $path) {
		if ($debug) {
			printf STDERR "Debug:  Called from \'%s\' at line \'%u\'...\n", $srcname, $subline;
		}
		printf STDERR "Error:  The path \'%s\' needs to be a directory!  Stopped\n", $path;
		return 0;
	}

	return 1;
}

##  Check that a filename has been given
sub CheckFile {
	my ($filename, $srcname, $subline, $debug) = @_;

	if (!-e $filename) {
		if (1) {
			printf STDERR "Debug:  Called from \'%s\' at line \'%u\'...\n", $srcname, $subline;
		}
		printf STDERR "Error:  The path \'%s\' does not exist!\n", $filename;
		return 0;
	}

	if (!-f $filename) {
		if ($debug) {
			printf STDERR "Debug:  Called from \'%s\' at line \'%u\'...\n", $srcname, $subline;
		}
		printf STDERR "Error:  The file \'%s\' does not exist!\n", $filename;
		return 0;
	}

	return 1;
}

##  Check that a program has been given
sub CheckProgram {
	my ($program, $srcname, $subline, $debug) = @_;

	if (!-e $program) {
		if (1) {
			printf STDERR "Debug:  Called from \'%s\' at line \'%u\'...\n", $srcname, $subline;
		}
		printf STDERR "Error:  The path \'%s\' does not exist!\n", $program;
		return 0;
	}

	if (!-x $program) {
		if ($debug) {
			printf STDERR "Debug:  Called from \'%s\' at line \'%u\'...\n", $srcname, $subline;
		}
		printf STDERR "Error:  The path \'%s\' is not an executable program!\n", $program;
		return 0;
	}

	return 1;
}


##############################################
#  Functions related to handling ORF reference
##############################################

sub CountORFs {
	my $orf_ref = shift;
	return scalar (@$orf_ref);
}


sub GetNameORFs {
	my ($pos, $orf_ref) = @_;

	my $orf_tmp = @$orf_ref[$pos];
	my @orf_array = split (/,/, $orf_tmp, 3);

	return $orf_array[0];
}


sub GetStartORFs {
	my ($pos, $orf_ref) = @_;

	my $orf_tmp = @$orf_ref[$pos];
	my @orf_array = split (/,/, $orf_tmp, 3);

	return $orf_array[1];
}


sub GetEndORFs {
	my ($pos, $orf_ref) = @_;

	my $orf_tmp = @$orf_ref[$pos];
	my @orf_array = split (/,/, $orf_tmp, 3);

	return $orf_array[2];
}

##  Returns 1 if there was an error in the ORF format
sub CheckORF {
	my $orf_tmp = shift;

	my ($orf_name, $orf_start, $orf_end, $empty) = split (/,/, $orf_tmp, 4);

	if (defined $empty) {
		printf STDERR "Error:  ORF format should be \'name,start,end\', with no commas in the name.\n";
		return 1;
	}

	if ($orf_start !~ /^\d+$/) {
		printf STDERR "Error:  ORF start position must be an integer.\n";
		return 1;
	}

	if ($orf_end !~ /^\d+$/) {
		printf STDERR "Error:  ORF end position must be an integer.\n";
		return 1;
	}

	return 0;
}


1;

