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

my $G_MAX_ATTEMPTS = 100;
my $G_TOLERANCE = 10;
my $G_OFFSET = 100;
my $G_PERCENT_VIRUS_LIMIT = 5;

my $fasta_width = 50;
my $chromosome_fn = "";
my $chromosome_seq = "";
my $chromosome_len = 0;
my $num_insert = 0;
my $freq_insert = 0;
my $virus_block_len = 0;

my $virus_fn = "";
my $virus_seq = "";
my $virus_len = 0;
my $low_virus_proportion = 0;
my $high_virus_proportion = 1.00;
my $low_virus_nts = 100;
my $high_virus_nts = 100;

my $is_gzip_file = "";
my $stats_insert_fn = "";
my $stats_fusion_fn = "";


##############################
##  Process arguments

##  Create AppConfig and AppConfig::Getopt objects
my $config = AppConfig -> new ({
        GLOBAL => {
            DEFAULT => undef,     ##  Default value for new variables
        }
    });

my $getopt = AppConfig::Getopt -> new ($config);

$config -> define ("chromosome", {
            ARGCOUNT => AppConfig::ARGCOUNT_ONE,
            ARGS => "=s",
        });                            ##  Filename to the chromosome
$config -> define ("num-insert", {
            ARGCOUNT => AppConfig::ARGCOUNT_ONE,
            ARGS => "=i",
        });                            ##  Number of insertions to perform
$config -> define ("freq-insert", {
            ARGCOUNT => AppConfig::ARGCOUNT_ONE,
            ARGS => "=i",
        });                            ##  Frequency of virus insertion
$config -> define ("virus-block-len", {
            ARGCOUNT => AppConfig::ARGCOUNT_ONE,
            ARGS => "=i",
        });                            ##  Length of a single virus "block"
$config -> define ("virus", {
            ARGCOUNT => AppConfig::ARGCOUNT_ONE,
            ARGS => "=s",
        });                            ##  Filename to the virus
$config -> define ("low-virus", {
            ARGCOUNT => AppConfig::ARGCOUNT_ONE,
            ARGS => "=i",
        });                            ##  Lower bound on percentage of virus
$config -> define ("high-virus", {
            ARGCOUNT => AppConfig::ARGCOUNT_ONE,
            ARGS => "=i",
        });                            ##  Upper bound on percentage of virus
$config -> define ("stats-insert", {
            ARGCOUNT => AppConfig::ARGCOUNT_ONE,
            DEFAULT => "insert.report",
            ARGS => "=s",
        });                            ##  Output file for statistics (insertion position)
$config -> define ("stats-fusion", {
            ARGCOUNT => AppConfig::ARGCOUNT_ONE,
            DEFAULT => "fusion.report",
            ARGS => "=s",
        });                            ##  Output file for statistics (exact location)
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

if (!defined $config -> get ("chromosome")) {
  printf STDERR "EE\tChromosome required with the --chromosome option.\n";
  exit (1);
}
$chromosome_fn = $config -> get ("chromosome");

if ((!defined $config -> get ("num-insert")) && (!defined $config -> get ("freq-insert")) && (!defined $config -> get ("virus-block-len"))) {
  printf STDERR "EE\tAt least one of --num-insert, --freq-insert, or --virus-block-len must be supplied.  No defaults are provided.\n";
  exit (1);
}

my $arg_count = 0;
if (defined $config -> get ("num-insert")) {
  $num_insert = $config -> get ("num-insert");
  $arg_count++;
}
elsif (defined $config -> get ("freq-insert")) {
  $freq_insert = $config -> get ("freq-insert");
  $arg_count++;
}
elsif (defined $config -> get ("virus-block-len")) {
  $virus_block_len = $config -> get ("virus-block-len");
  $arg_count++;
}

if ($arg_count != 1) {
  printf STDERR "EE\tOnly one of --num-insert, --freq-insert, or --virus-block-len can be supplied.\n";
  exit (1);
}

if (!defined $config -> get ("virus")) {
  printf STDERR "EE\tVirus required with the --virus option.\n";
  exit (1);
}
$virus_fn = $config -> get ("virus");

##  Check the lower and upper bounds of the virus proportion
##    Four cases:
##      1)  Neither lower and upper not set --> Error
##      2)  Both lower and upper are set
##          a)  Both lower and upper have the same value
##          b)  Lower and upper have different values
##      3)  Only lower is set -- Set higher to the same value
##      4)  Only upper is set -- Set lower to the same value
if ((!defined $config -> get ("low-virus")) && (!defined $config -> get ("high-virus"))) {
  printf STDERR "EE\tAt least one of --low-virus or --high-virus must be supplied.  No defaults are provided.\n";
  exit (1);
}
elsif ((defined $config -> get ("low-virus")) && (defined $config -> get ("high-virus"))) {
  if (($config -> get ("low-virus")) > ($config -> get ("high-virus"))) {
    printf STDERR "EE\tThe value for --low-virus must be smaller than the value for --high-virus.\n";
    exit (1);
  }
  $low_virus_proportion = $config -> get ("low-virus");
  $high_virus_proportion = $config -> get ("high-virus");
}
elsif (defined $config -> get ("low-virus")) {
  $low_virus_proportion = $config -> get ("low-virus");
  $high_virus_proportion = $config -> get ("low-virus");
}
elsif (defined $config -> get ("high-virus")) {
  $low_virus_proportion = $config -> get ("high-virus");
  $high_virus_proportion = $config -> get ("high-virus");
}

##  At this point, both $low_virus_proportion and $high_virus_proportion have been set but are still percentages
if (($low_virus_proportion <= 0) || ($high_virus_proportion > 100)) {
  printf STDERR "EE\tThe range bounded by --low-virus and --high-virus must be in (0, 100].\n";
  exit (1);
}

##  Convert the percentages to proportions
$low_virus_proportion = $config -> get ("low-virus") / 100;
$high_virus_proportion = $config -> get ("high-virus") / 100;

if (defined $config -> get ("gzip")) {
  $is_gzip_file = 1;
}

if (!defined $config -> get ("stats-insert")) {
  printf STDERR "EE\tFilename for the insertion statistics needed with the --stats-insert parameter!\n";
  exit (-1);
}
$stats_insert_fn = $config -> get ("stats-insert");

if (!defined $config -> get ("stats-fusion")) {
  printf STDERR "EE\tFilename for the location of the fusion position statistics needed with the --stats-fusion parameter!\n";
  exit (-1);
}
$stats_fusion_fn = $config -> get ("stats-fusion");


##############################
##  Print information out

printf STDERR "II\tRandom seed to be used:  %u\n", srand ();


##############################
##  Read the virus in

my $fp = OpenFastqFileRead ($virus_fn, $is_gzip_file);
my $virus_header = GetNextFastqLine ($fp, $is_gzip_file);
while (1) {
  my $line = GetNextFastqLine ($fp, $is_gzip_file);
  if (length ($line) == 0) {
    last;
  }
  chomp $line;

  ##  If there is more than one FASTA sequence, just take the first one
  if ($line =~ /\>/) {
    last;
  }

  $virus_seq = $virus_seq.$line;
}
$virus_len = length ($virus_seq);

CloseFastqFile ($fp, $is_gzip_file);

##  Handle variables related to virus' length or virus' block length
my $tmp_block_length = $virus_len;
if (defined $config -> get ("virus-block-len")) {
  $tmp_block_length = $virus_block_len;
}

if ($config -> get ("low-virus") == 100) {
  $low_virus_nts = $tmp_block_length;
}
else {
  $low_virus_nts = int ($tmp_block_length * $low_virus_proportion);
}

if ($config -> get ("high-virus") == 100) {
  $high_virus_nts = $tmp_block_length;
}
else {
  $high_virus_nts = int ($tmp_block_length * $high_virus_proportion);
}


if ($config -> get ("verbose")) {
  if (defined $config -> get ("virus-block-len")) {
    printf STDERR "II\tLength of a (non-redundant) virus block:  %u\n", $virus_block_len;
    printf STDERR "II\t  Number of insertions:\t%u\n", int ($virus_len / $virus_block_len);
  }
  printf STDERR "II\tVirus sequence length:\t%u\n", $virus_len;
  printf STDERR "II\tNucleotides from virus:\t[%u, %u]\n", $low_virus_nts, $high_virus_nts;
  if ($low_virus_nts == 0) {
    printf STDERR "II\tThe low virus length cannot be 0.\n";
  }
}


##############################
##  Read the chromosome in

$fp = OpenFastqFileRead ($chromosome_fn, $is_gzip_file);
my $chromosome_header = GetNextFastqLine ($fp, $is_gzip_file);
if ($chromosome_header =~ /^>(.+)$/) {
  $chromosome_header = $1;
}

while (1) {
  my $line = GetNextFastqLine ($fp, $is_gzip_file);
  if (length ($line) == 0) {
    last;
  }
  chomp $line;

  ##  If there is more than one FASTA sequence, just take the first one
  if ($line =~ /\>/) {
    last;
  }

  $chromosome_seq = $chromosome_seq.$line;
}
$chromosome_len = length ($chromosome_seq);

CloseFastqFile ($fp, $is_gzip_file);


##  Handle variables related to chromosome's length
if (defined $config -> get ("virus-block-len")) {
  $num_insert = int ($virus_len / $virus_block_len)
}
elsif ($num_insert == 0) {
  $num_insert = int ($chromosome_len / $freq_insert);
}

if ($config -> get ("verbose")) {
  printf STDERR "II\tChromosome sequence length:\t%u\n", $chromosome_len;
  printf STDERR "II\tNumber of viruses to insert:\t%u\n", $num_insert;
}
if ($num_insert > $chromosome_len * $G_PERCENT_VIRUS_LIMIT) {
  printf STDERR "EE\tThere can be no more than %.1f%% of the chromosome with virus insertions\nEE\t  (i.e., no more than %u).\nEE\t  Please reduce --num-insert or increase the value to --freq-insert!\n", $G_PERCENT_VIRUS_LIMIT, int ($chromosome_len * $G_PERCENT_VIRUS_LIMIT);
  exit (1);
}


##############################
##  Decide positions for virus insertions

if ($config -> get ("verbose")) {
  printf STDERR "II\tInitializing arrays...\n";
}

##  Waste of space, but used to indicate places to insert viruses;
##    Basically, an array of binary flags
my @chromosome_mark;

##  Used to indicate the virus length added at each position in the chromosome
##    Also serves as a binary flag.
my @chromosome_virus_len;

##  Initialize
for (my $i = 0; $i < $chromosome_len; $i++) {
  $chromosome_mark[$i] = 0;
  $chromosome_virus_len[$i] = 0;
}

if ($config -> get ("verbose")) {
  printf STDERR "II\tDetermining where to insert viruses...\n";
}

my $j = 0;
my $attempts = 0;
while ($j < $num_insert) {
  my $pos = int (rand ($chromosome_len));

  if ($attempts >= $G_MAX_ATTEMPTS) {
    printf STDERR "EE\tThere have been too many failed attempts in finding an insertion site...\n";
    exit (1);
  }

  ##  If the position is close to the ends of the chromosome, search again
  if (($pos - $G_TOLERANCE < 0) || ($pos + $G_TOLERANCE > $chromosome_len)) {
    $attempts++;
    next;
  }

  ##  If the area around the insertion point has at least one N, then try again
  my $read_seq  = substr ($chromosome_seq, $pos - $G_TOLERANCE, 2 * $G_TOLERANCE);
  if ($read_seq =~ /N/) {
    $attempts++;
    next;
  }

  $attempts = 0;
  if ($chromosome_mark[$pos] == 0) {
    $chromosome_mark[$pos] = 1;
    $j++;
  }
}


##############################
##  Perform the actual insertion

if ($config -> get ("verbose")) {
  printf STDERR "II\tPerforming viral insertion...\n";
}

##  Open output files for the statistics
open (my $stats_insert_fp, ">", $stats_insert_fn) or die "EE\tCould not open $stats_insert_fn for output.\n";
open (my $stats_fusion_fp, ">", $stats_fusion_fn) or die "EE\tCould not open $stats_fusion_fn for output.\n";

my @blocks;
if (defined $config -> get ("virus-block-len")) {
  for (my $i = 0; $i < $num_insert; $i++) {
    $blocks[$i] = 0;
  }
}

my $block_id = 0;
my $total_virus_len = 0;
for (my $i = $chromosome_len - 1; $i >= 0; $i--) {
  if ($chromosome_mark[$i] == 1) {
    my $pos = $i;

    ##  Determine the virus to insert:  (1)  Length and (2) Position in the virus.
    my $sub_virus_len = int (rand ($high_virus_nts - $low_virus_nts)) + $low_virus_nts;
    my $sub_virus_start = 0;
    ##  Starting position of the virus is either:  (1) somewhere within a random block or
    ##    (2) somewhere within the entire virus.
    if (defined $config -> get ("virus-block-len")) {
      my $block_id = 0;
      do {
        $block_id = int (rand ($num_insert));
        $attempts++;
        if ($attempts == $G_MAX_ATTEMPTS) {
          printf STDERR "EE\tThere have been too many failed attempts in finding a starting viral position...\n";
          exit (1);
        }
      } while ($blocks[$block_id] == 1);
      $attempts = 0;
      $blocks[$block_id] = 1;

      $sub_virus_start = int (rand ($virus_block_len - $sub_virus_len)) + ($block_id * $virus_block_len);
    }
    else {
      $sub_virus_start = int (rand ($virus_len - $sub_virus_len));
    }
    my $sub_virus_end = $sub_virus_start + $sub_virus_len;
    my $sub_virus_seq = substr ($virus_seq, $sub_virus_start, $sub_virus_len);

    ##  Insert virus
    substr ($chromosome_seq, $pos, 0) = $sub_virus_seq;

    ##  Report the results
    printf $stats_insert_fp "%s\t", $chromosome_header;
    printf $stats_insert_fp "%u\t", $pos;  ##  Position in the chromosome that this virus was inserted into
    printf $stats_insert_fp "%u\t", $sub_virus_start;
    printf $stats_insert_fp "%u\t", $sub_virus_end;
    printf $stats_insert_fp "%u\n", $sub_virus_len;

    ##  Record the virus length inserted at each position
    $chromosome_virus_len[$pos] = $sub_virus_len;

    ##  Keep track of how much of the virus we inserted
    $total_virus_len += $sub_virus_len;
  }
}

##  Get the length of the new chromosome
my $new_chromosome_len = length ($chromosome_seq);

##  Print the fusion statistics
my $viruses_sofar = 0;
for (my $i = 0; $i < $chromosome_len; $i++) {
  if ($chromosome_virus_len[$i] != 0) {
    ##  Print the left junction
    my $offset = $G_OFFSET;
    my $midpoint = $i + $viruses_sofar;
    if ($offset > $midpoint) {
      $offset = $midpoint;
    }
    my $pre_fusion = substr ($chromosome_seq, $midpoint - $offset, $offset);
    my $post_fusion = substr ($chromosome_seq, $midpoint, $offset);

    printf $stats_fusion_fp "%s\t", $chromosome_header;
    printf $stats_fusion_fp "%u\t", $midpoint;
    printf $stats_fusion_fp "L\t";
    printf $stats_fusion_fp "%s\t%s\n", $pre_fusion, $post_fusion;

    ##  Print the right junction
    $offset = $G_OFFSET;
    $midpoint = $i + $viruses_sofar + $chromosome_virus_len[$i];
    if ($offset > ($new_chromosome_len - $midpoint)) {
      ##  Ensure the flanking region does not go off the end of the chromosome
      $offset = $new_chromosome_len - $midpoint;
    }

    $pre_fusion = substr ($chromosome_seq, $midpoint - $offset, $offset);
    $post_fusion = substr ($chromosome_seq, $midpoint, $offset);

    printf $stats_fusion_fp "%s\t", $chromosome_header;
    printf $stats_fusion_fp "%u\t", $midpoint;
    printf $stats_fusion_fp "R\t";
    printf $stats_fusion_fp "%s\t%s\n", $pre_fusion, $post_fusion;

    $viruses_sofar += $chromosome_virus_len[$i];
  }
}


##  Close the output files for the statistics
close ($stats_insert_fp);
close ($stats_fusion_fp);

##  Add new lines to the chromosome
$chromosome_seq =~ s/(.{$fasta_width})/$1\n/g;

if ($config -> get ("gzip")) {
  my $out_fh = gzopen (\*STDOUT, "wb") or die "Can't open stdout for writing : $gzerrno";
  $chromosome_header = ">".$chromosome_header."\n";
  $out_fh -> gzwrite ($chromosome_header);
  $out_fh -> gzwrite ($chromosome_seq);
  $out_fh -> gzwrite ("\n");
  $out_fh -> gzclose ();
}
else {
  printf ">%s\n", $chromosome_header;
  printf "%s\n", $chromosome_seq;
}

if ($config -> get ("verbose")) {
  printf STDERR "II\tOld chromosome sequence length:\t%u\n", $chromosome_len;
  printf STDERR "II\tNew chromosome sequence length:\t%u\n", $new_chromosome_len;
  printf STDERR "II\tTotal virus length:\t%u\n", $total_virus_len;
  printf STDERR "II\tPercentage of virus:\t%.3f%%\n", $total_virus_len / $new_chromosome_len * 100;
  printf STDERR "\n";
  printf STDERR "II\tReports have been sent to:\n";
  printf STDERR "II\t  Insertions:  %s\n", $stats_insert_fn;
  printf STDERR "II\t  Fusions:  %s\n", $stats_fusion_fn;
  printf STDERR "II\tSee this script's \'perldoc\' for a description of the fields.\n";
}


=pod

=head1 NAME

simulate-viralfusion.pl -- Simulate viral fusions.

=head1 SYNOPSIS

B<simulate-viralfusion.pl> [OPTIONS]

=head1 DESCRIPTION

Given a viral genome and a chromosome, portions of the viral genome is inserted into the chromosome at random positions.  The part of the viral genome that is added is also decided at random.  No mutation of either is performed.  At most 5% of the chromosome can have viral insertions; to change this, change the value to $G_PERCENT_VIRUS_LIMIT.  When an insertion point in the chromosome is selected, $G_TOLERANCE bases are checked up and down stream to make sure there are no N's.  If there is at least one N in either region, then a new insertion point is selected.

The new chromosome sequence is sent to STDOUT.  Information about where the insertions were made are sent to two different files:  an insertion report and a fusion report.  Both are in tab-separated format with 0-based coordinates.

The insertion report indicates where in the chromosome insertions occurred; this is using co-ordinates relative to the original chromosome.  Because of this, the file itself is sorted in reverse due to how the insertions were added (from the end instead of from the beginning).

The insertion report consists of the following 5 fields:

=over 5

=item 1. chromosome name

=item 2. position of the insertion into the chromosome

=item 3. start position of the virus (inclusive)

=item 4. end position of the virus (exclusive)

=item 5. length of the part of the virus that was inserted -- equals to (field #4 - field #3)

=back

The fusion report indicates where the fusion junctions actually occurred using the coordinate system of the infected genome.  This is the two end-points of each viral insertion.  For n viral insertions, there are 2n fusion junctions.

The fusion report has 3 fields:

=over 5

=item 1. chromosome name

=item 2. left position of the viral insertion OR one position to the right of the viral insertion

=item 3. "L" or "R", depending on left or right position

=item 4. left flanking region around the junction site; length is -$G_OFFSET base pairs

=item 5. right flanking region around the junction site; length is +$G_OFFSET base pairs

=back

As an example, the following indicates virus (digits) being inserted into a chromosome of A's.  The left and right positions are indicated below with `L' and `R', respectively.

=over 10

>sample chromosome

2345AAAAAAAA5678AAAA45678AAAAAAAA

0         1         2         3
012345678901234567890123456789012
L   R       L   R   L    R

=back

Note that at the boundary between virus-human or human-virus, the position of the junction is indicated using the position immediately to the B<right> of the junction.  In the example above, the junction sites are 0, 4, 12, 16, 20, and 25.

=head1 OPTIONS

=over 5

=item --chromosome F<path>

Path to the chromosome, in FASTA format.  Only the first FASTA sequence in this file is taken.

=item --num-insert I<integer>

Number of virus insertions to perform.  Only one of --num-insert, --freq-insert, or --virus-block-len is allowed.

=item --freq-insert I<integer>

Frequency of virus insertions.  The actual number that is inserted is determined by (length of chromosome / frequency of insertions).  Only one of --num-insert, --freq-insert, or --virus-block-len is allowed.

=item --virus-block-len I<integer>

Divide the virus into non-overlapping blocks of this size.  The inserted viruses are randomly selected from each of these blocks.  Thus, the number of insertions is equal to the number of blocks that can be created.  Only one of --num-insert, --freq-insert, or --virus-block-len is allowed.

=item --virus F<path>

Path to the virus genome, in FASTA format.  Only the first FASTA sequence in this file is taken.

=item --low-virus I<integer>

Lower bound on the B<percent> of the virus to take.  If you do not want to provide a range, then give either --low-virus or --high-virus.

=item --high-virus I<integer>

Upper bound on the B<percent> of the virus to take.  Set this to 100 if you always want to take all of the virus.  If you do not want to provide a range, then give either --low-virus or --high-virus.

=item --stats-insert I<file>

File where the insertion report is written to.  Default is F<insert.report>.

=item --stats-fusion I<file>

File where the fusion junction report is written to.  Default is F<fusion.report>.

=item --gzip

Employ gzip for both input and output files.  That is, the FASTA genomes and the infected FASTA chromosome, respectively.

=item --verbose

Display verbose information about the execution of the program.

=item --help

Display this help message.

=back

=head1 EXAMPLES

The following will insert the HIV virus into the H. sapiens' chromosome 1.  25-50% of the virus is inserted each time.  The number of insertions is determined by (length of chromosome 1 / 1000000).  This does B<not> mean that the virus is inserted at every 1,000,000 nucleotides.  Input/output is done using compression with gzip.  Verbose output is requested by the user.

=over 5

./simulate-viralfusion.pl --chromosome chr1.fa.gz --virus NC_001802.1.fa.gz --gzip --freq-insert 1000000 --low-virus 25 --high-virus 50 --verbose 1>new_chr1.fa

=back

The following will insert the HIV virus into the H. sapiens' chromosome 1.  The virus is segmented into non-overlapping "blocks" such that each one is sampled B<once> independently.  75-100% of each block is inserted each time (or 300-400 nts in this example).  The number of insertions is determined by (length of viral sequence / virus block length).  Input/output is done using compression with gzip.  Verbose output is requested by the user.

=over 5

./simulate-viralfusion.pl --chromosome chr1.fa.gz --virus NC_001802.1.fa.gz --gzip --virus-block-len 400 --low-virus 75 --high-virus 100 --verbose 1>new_chr1.fa.gz

=back

=head1 AUTHOR

Raymond Wan <rwan@cuhk.edu.hk>

=head1 COPYRIGHT

Copyright (C) 2012, Raymond Wan, All rights reserved.

