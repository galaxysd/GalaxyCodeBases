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

my $chromosome_size = 0;
my $num_reads = 0;
my $read_len = 0;
my $coverage = 0;

my $num_to_select = 0;
my $num_selected = 0;
my @reads_selected;


##############################
##  Process arguments

##  Create AppConfig and AppConfig::Getopt objects
my $config = AppConfig -> new ({
        GLOBAL => {
            DEFAULT => undef,     ##  Default value for new variables
        }
    });

my $getopt = AppConfig::Getopt -> new ($config);

$config -> define ("chrsize", {
            ARGCOUNT => AppConfig::ARGCOUNT_ONE,
            ARGS => "=i",
        });                            ##  Chromosome size
$config -> define ("readlen", {
            ARGCOUNT => AppConfig::ARGCOUNT_ONE,
            ARGS => "=i",
        });                            ##  Read length
$config -> define ("numreads", {
            ARGCOUNT => AppConfig::ARGCOUNT_ONE,
            ARGS => "=i",
        });                            ##  Number of reads
$config -> define ("coverage", {
            ARGCOUNT => AppConfig::ARGCOUNT_ONE,
            ARGS => "=i",
        });                            ##  Desired coverage
$config -> define ("single!", {
            DEFAULT => 0,
        });                            ##  Single-ended reads
$config -> define ("gzip!", {
            DEFAULT => 0,
        });                            ##  Employ gzip
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

if (!defined $config -> get ("chrsize")) {
  printf STDERR "EE\tChromosome size required with the --chrsize option.\n";
  exit (1);
}
$chromosome_size = $config -> get ("chrsize");

if (!defined $config -> get ("readlen")) {
  printf STDERR "EE\tRead length required with the --readlen option.\n";
  exit (1);
}
$read_len = $config -> get ("readlen");

if (!defined $config -> get ("numreads")) {
  printf STDERR "EE\tNumber of reads required with the --numreads option.\n";
  exit (1);
}
$num_reads = $config -> get ("numreads");

if (!defined $config -> get ("coverage")) {
  printf STDERR "EE\tCoverage required with the --coverage option.\n";
  exit (1);
}
$coverage = $config -> get ("coverage");


##############################
##  Determine number of reads that need to be selected
$num_to_select = ($coverage * $chromosome_size);
if ($config -> get ("single")) {
  $num_to_select = int ($num_to_select / ($read_len));
}
else {
  $num_to_select = int ($num_to_select / (2 * $read_len));
}

if ($config -> get ("verbose")) {
  printf STDERR "II\tNumber of reads to select:  %u\n", $num_to_select;
  printf STDERR "II\tNumber of reads available:  %u\n", $num_reads;
}

if ($num_to_select > $num_reads) {
  printf STDERR "EE\tNot enough reads available!\n";
  exit (1);
}


##############################
##  Initialization and select reads

##  Use <= since we want $num_reads to be included
for (my $i = 0; $i <= $num_reads; $i++) {
  $reads_selected[$i] = 0;
}

for (my $i = 0; $i < $num_to_select; $i++) {
  my $pos = 0;
  do {
    $pos = rand ($num_reads) + 1;  ##  Need to add 1 since the search is 0-based.
  } while ($reads_selected[$pos] == 1);
  $reads_selected[$pos] = 1;
}


##############################
##  Send numbers to compressed gzip file

if ($config -> get ("gzip")) {
  my $fh = gzopen (\*STDOUT, "wb") or die "Can't open STDOUT for reading : $gzerrno";
  for (my $i = 0; $i <= $num_reads; $i++) {
    if ($reads_selected[$i] == 1) {
      my $outline = $i."\n";
      $fh -> gzwrite ($outline);
    }
  }
  $fh -> gzclose ();
}
else {
  for (my $i = 0; $i <= $num_reads; $i++) {
    if ($reads_selected[$i] == 1) {
      printf "%u\n", $i;
    }
  }
}


=pod

=head1 NAME

sample-reads.pl -- "Sample" reads by output read IDs.

=head1 SYNOPSIS

B<sample-reads.pl> [OPTIONS]

=head1 DESCRIPTION

This script randomly "samples" reads using a uniform distribution.  Actually reads are not read by this script.  Instead, it prints a list of numbers to STDOUT which is used by other scripts for selecting reads.  How many reads are selected depends on:

=over 5

=item 1. size of the chromosome,

=item 2. read length,

=item 3. coverage,

=item 4. number of available reads, and

=item 5. whether single or paired-end reads are required.

=back

This is because of this relationship:

=over 5

number of reads = (coverage * chromosome size) / (read length)

=back

(The read length is multiplied by 2 for paired-end reads.  An error message is returned if the number of reads needed is less than the number available.

=head1 OPTIONS

=over 5

=item --chrsize I<integer>

Size of the chromosome.

=item --chrsize I<readlen>

The length of the reads desired.

=item --chrsize I<numreads>

The actual number of reads available.

=item --chrsize I<coverage>

The coverage desired.

=item --single

Single-ended reads are required.  The default is paired-end reads.

=item --gzip

Employ the gzip library for compressing the output.

=item --verbose

Display verbose information about the execution of the program.

=item --help

Display this help message.

=back

=head1 EXAMPLES

The following will generate a list of numerical read IDs assuming a chromosome size of 249,250,621 nucleotides (UCSC hg19 chromosome 1), a read length of 101 base pairs, and a coverage of 100x.  There are 250,000,000 reads to select from.  The output is sent as text to STDOUT, with one value per line, in sorted order.

=over 5

./sample-reads.pl --chrsize 249250621 --readlen 101 --numreads 250000000 --coverage 100

=back

=head1 AUTHOR

Raymond Wan <rwan@cuhk.edu.hk>

=head1 COPYRIGHT

Copyright (C) 2012, Raymond Wan, All rights reserved.


