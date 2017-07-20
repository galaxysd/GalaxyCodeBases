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

##  $LastChangedDate: 2012-11-27 17:32:45 +0800 (Tue, 27 Nov 2012) $
##  $LastChangedRevision: 1092 $

use diagnostics;
use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/../include";

##  Module for general subroutines; functions to be used have "Fastq" in the name,
##  but the ones being used are not specific to FASTQ; FASTA is ok, too.
use General; 

##  Modules not provided with this program
use AppConfig;
use AppConfig::Getopt;
use Pod::Usage;
use Switch;
use Cwd;
use Compress::Zlib;


##############################
##  Global variables

my $param = 0;  ##  Parameter used for printing; can be any integer
my $fusion_report_fn = "";
my $cs_results_fn = "";
my $cs_window_size = 0;
my $rp_results_fn = "";

my %cs_results;
my %rp_results;

##  Hash that maps based on position in the chromosome.
##    0 -- no meaning; never used
##    1 -- means that this is a junction position
##    >1 -- junction position that's been found
my %fusion_pos_hash;


##############################
##  Process arguments

##  Create AppConfig and AppConfig::Getopt objects
my $config = AppConfig -> new ({
        GLOBAL => {
            DEFAULT => undef,     ##  Default value for new variables
        }
    });

my $getopt = AppConfig::Getopt -> new ($config);

$config -> define ("param", {
            ARGCOUNT => AppConfig::ARGCOUNT_ONE,
            ARGS => "=i",
            DEFAULT => 0,
        });                            ##  Parameter used for printing
$config -> define ("fusion-report", {
            ARGCOUNT => AppConfig::ARGCOUNT_ONE,
            ARGS => "=s",
        });                            ##  Filename of the fusion report
$config -> define ("cs-results", {
            ARGCOUNT => AppConfig::ARGCOUNT_ONE,
            ARGS => "=s",
        });                            ##  Filename of the CS results
$config -> define ("cs-window-size", {
            ARGCOUNT => AppConfig::ARGCOUNT_ONE,
            ARGS => "=i",
            DEFAULT => 0,
        });                            ##  Window size for the CS method
$config -> define ("rp-results", {
            ARGCOUNT => AppConfig::ARGCOUNT_ONE,
            ARGS => "=s",
        });                            ##  Filename of the RP results
$config -> define ("verbose!", {
            DEFAULT => 0,
        });                            ##  Verbose output
$config -> define ("help!", {
            DEFAULT => 0,
        });                            ##  Help screen

##  Process the command-line options
$config -> getopt ();


##############################
##  Validate the settings
if ($config -> get ("help")) {
  pod2usage (-verbose => 0);
  exit (1);
}


$param = $config -> get ("param");


if (!defined $config -> get ("fusion-report")) {
  printf STDERR "EE\tFilename to the fusion report.\n";
  exit (1);
}
$fusion_report_fn = $config -> get ("fusion-report");


if ((!defined $config -> get ("cs-results")) && (!defined $config -> get ("rp-results"))) {
  printf STDERR "EE\tEither or both --cs-results and --rep-results is required.\n";
  exit (1);
}
if (defined $config -> get ("cs-results")) {
  $cs_results_fn = $config -> get ("cs-results");
  if (!defined $config -> get ("cs-window-size")) {
    printf STDERR "EE\tWindow for the CS method required with the --cs-window-size option.\n";
    exit (1);
  }
  $cs_window_size = $config -> get ("cs-window-size");
}
if (defined $config -> get ("rp-results")) {
  $rp_results_fn = $config -> get ("rp-results");
}


##############################
##  Initialization followed by open fusion report

open (my $fp, "<", $fusion_report_fn) or die "EE\tCould not open $fusion_report_fn for reading.\n";
while (<$fp>) {
  my $line = $_;
  chomp $line;

  my ($chr, $pos, $type, $left_seq, $right_seq) = split /\t/, $line;
  if (length ($left_seq) != length ($right_seq)) {
    printf STDERR "EE\tError in reading in fusion report.  Left and right sequences do not match in length.\n";
    exit (1);
  }

  $fusion_pos_hash{$pos} = 1;

#   printf STDERR "DD\tSet position %u from fusion report...\n", $pos;
}
close ($fp);


##############################
##  Read in CS method
if (defined $config -> get ("cs-results")) {
  if ($config -> get ("verbose")) {
    printf STDERR "II\tProcessing CS results.\n";
  }

  open (my $fp, "<", $cs_results_fn) or die "EE\tCould not open $cs_results_fn.\n";
  while (<$fp>) {
    my $line = $_;
    chomp $line;

    my ($read_id, $pos) = split /\t/, $line;
    if (($read_id !~ /^\d+$/) || ($pos !~ /^\d+$/)) {
      printf STDERR "EE\tProblem parsing \'%s\' in CS results file.\n", $line;
      exit (1);
    }

    $cs_results{$read_id} = $pos;
#     printf STDERR "DD\tAssigned read ID [%u] --> %u\n", $read_id, $pos;
  }
  close ($fp);
}


##############################
##  Read in RP method and set array
if (defined $config -> get ("rp-results")) {
  if ($config -> get ("verbose")) {
    printf STDERR "II\tProcessing RP results.\n";
  }

  open (my $fp, "<", $rp_results_fn) or die "EE\tCould not open $rp_results_fn.\n";
  while (<$fp>) {
    my $line = $_;
    chomp $line;

    my ($read_id, $start, $end) = split /\t/, $line;
    if (($read_id !~ /^\d+$/) || ($start !~ /^\d+$/) || ($end !~ /^\d+$/)) {
      printf STDERR "EE\tProblem parsing \'%s\' in RP results file.\n", $line;
      exit (1);
    }

    $rp_results{$read_id} = $start."\t".$end;
  }
  close ($fp);
}


##############################
##  Read in read list in plain text
if ($config -> get ("verbose")) {
  printf STDERR "II\tReading list of reads.\n";
}


my $cs_found = 0;
my $rp_found = 0;
my $inaccuracy = 0;
while (<STDIN>) {
  my $read_id = $_;
  chomp $read_id;

  if ($read_id !~ /^\d+$/) {
    printf STDERR "EE\tInvalid value \'%s\' in the list of reads IDs.\n", $read_id;
    exit (1);
  }

#   printf STDERR "DD\tRead ID [%u] read.\n", $read_id;

  ##  Translate read ID into a position, if available
  if (defined $cs_results{$read_id}) {
    $cs_found = 0;
    my $found_pos = 0;
    my $cs_pos = $cs_results{$read_id};

#     printf STDERR "DD\tRead ID %u corresponds to position %u\n", $read_id, $cs_pos;

    ##  Perfect hit
    if ((defined $fusion_pos_hash{$cs_pos}) && ($fusion_pos_hash{$cs_pos} >= 1)) {
      $found_pos = $cs_pos;
    }
    elsif ($cs_window_size != 0) {
      ##  In the case of CS results, search the area around the junction position
      my $end = $cs_pos + $cs_window_size;
      ##  Set the positions ahead
      for (my $i = $cs_pos; $i < $end; $i++) {
        if (defined $fusion_pos_hash{$i}) {
          $found_pos = $i;
          last;
        }
      }

      ##  Still haven't found it, so search backwards
      if ($found_pos == 0) {
        my $start = $cs_pos - $cs_window_size;
        if ($start < 0) {
          $start = 0;
        }
        ##  Set the positions behind
        for (my $i = $start; $i < $cs_pos; $i++) {
          if (defined $fusion_pos_hash{$i}) {
            $found_pos = $i;
            last;
          }
        }
      }
    }

    if ($found_pos != 0) {
      $fusion_pos_hash{$found_pos}++;
      $cs_found = 1;
    }
    else {
      $cs_found = 0;
    }
  }

  ##  Translate read ID into a start and end positions, if available
  if (defined $rp_results{$read_id}) {
    $rp_found = 0;
    my ($rp_start, $rp_end) = split /\t/, $rp_results{$read_id};

    ##  We search forward and reverse; if they meet at different positions, then we have a problem.
    for (my $i = $rp_start; $i < $rp_end; $i++) {
      if (defined $fusion_pos_hash{$i}) {
        $rp_found = 1;
        $fusion_pos_hash{$i}++;
        ##  We do not 'break' out of here because within the gap, there will be more
        ##  than two fusion junctions.
      }
    }
  }

  if (($cs_found == 0) && ($rp_found == 0)) {
    $inaccuracy++;
  }
}


##############################
##  Output percentage

my $accuracy = 0;
my $total = 0;
foreach my $key (sort (keys %fusion_pos_hash)) {
  if ($fusion_pos_hash{$key} > 1) {
    $accuracy++;
  }
  $total++;
}
printf STDOUT "XXX\t%u\t%u\t%u\t%u\n", $config -> get ("param"), $accuracy, $inaccuracy, $total;


=pod

=head1 NAME

evaluate-viralfusion.pl -- Evaluate the output from VFS.

=head1 SYNOPSIS

B<evaluate-viralfusion.pl> [OPTIONS]

=head1 DESCRIPTION

This script assists in the evaluation of the output from the VFS pipeline.  It combines the information about the true locations of the junction locations with the outputs of the CS or RP modules (either or both is required).  It uses a list of reads to yield an accuracy -- i.e., how many of the junction sites can be predicted with the given set of reads.  Read IDs are assumed to be I<numerical>.

In short, it does the following:

=over 5

=item 1. Reads in all of the fusion junction locations, whose format is described in the perldoc of F<simulate-viralfusion.pl>.

=item 2. Reads in the output from the CS module (if provided), which is a tab separated list of read IDs and locations on the infected genome.

=item 3. Reads in the output of the RP module (if provided), which is also tab separated and has the read ID, end position of the forward read, and start position of the reverse read (thus, the junction is in between these two positions).

=item 4. Reads in the list of reads.

=item 5. Determines if each read represents a viral junction using first the CS module and then the RP module.

=back

There is no harm if both modules represents the same (or different!) junction.  The final result reported is the number of junctions accounted for; not how many times each was accounted for.  The reads identified by the RP module cover a range of base pairs (the insertion size); the CS module can take a parameter to cover a window of bases from the specified position.  Thus, there is a chance more than one junction is identified by a read.

The evaluation process does B<not> aggregate reads to get a consensus.  Each read is processed independently of all other reads.

=head1 OPTIONS

=over 5

=item --param I<integer>

A numerical parameter that is printed with the output.  Used only for display and not used in any calculation.

=item --fusion-report F<file>

Path to the fusion report.

=item --cs-results F<file>

Path to the CS results.

=item --cs-window-size I<integer>

The window size to use with the data from the CS module.  Values are +/- from the predicted position.  Default is 0.

=item --rp-results F<file>

Path to the RP results.

=item --verbose

Display verbose information about the execution of the program.

=item --help

Display this help message.

=back

=head1 EXAMPLES

Evaluate the accuracy of the reads in the sample.txt file using the provided fusion report and the CS results with a window size of 0 (the default).

=over 5

cat sample.txt | ./evaluate-viralfusion.pl --fusion-report fusion.report --cs-results cs-results.txt --verbose >output.txt

=back

=head1 AUTHOR

Raymond Wan <rwan@cuhk.edu.hk>

=head1 COPYRIGHT

Copyright (C) 2012, Raymond Wan, All rights reserved.


