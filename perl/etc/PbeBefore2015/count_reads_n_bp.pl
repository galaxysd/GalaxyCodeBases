#!/usr/bin/perl
use strict;
use warnings;
use threads;

my %file;
while (<>) {
  chomp;
  $file{$_} = openfile($_);
}
print STDERR "Read file list complete!\n";

my %thread;
foreach (sort keys %file) {
  if (/sam.gz$/ or /bam$/) {
    $thread{$_} = threads->new(\&countsam, $file{$_});
  } elsif (/fq.gz$/ or /addN.gz$/) {
    $thread{$_} = threads->new(\&countfq, $file{$_});
  }
}

my %count;
foreach (sort keys %thread) {
  $count{$_} = $thread{$_}->join;
  print STDERR "Count $_ complete!\n";
}

foreach (sort keys %count) {
  print $_, "\t", (join "\t", @{$count{$_}}), "\n";
}

sub openfile {
  my $filename = shift;
  my $infile;
  if ($filename =~ /.bam$/) {
    open $infile, "-|", "samtools view $filename" or die "Error opening $filename: $!\n";
  } elsif ($filename =~ /.sam.gz$/) {
    open $infile, "-|", "samtools view -S $filename" or die "Error opening $filename: $!\n";
  } else {
    open $infile, "-|", "gzip -dc $filename" or die "Error opening $filename: $!\n";
  }
  return $infile;
}

sub countfq {
  my $in = shift;
  my ($reads, $bp) = (0, 0);
  while (<$in>) {
    chomp (my $seq = <$in>);
    <$in>;
    <$in>;
    ++$reads;
    next if $seq eq "N";
    my $length = length $seq;
    $bp += $length;
  }
  close $in;
  return [$reads, $bp];
}

sub countsam {
  my $in = shift;
  my ($aligned_reads, $aligned_bp) = (0, 0);
  while (<$in>) {
    my @a = split /\t/;
    ++$aligned_reads if $a[3];
    my @b = $a[5] =~ /(\d+)[MI]/g;
    $aligned_bp += $_ foreach @b;
  }
  close $in;
  return [$aligned_reads, $aligned_bp];
}
