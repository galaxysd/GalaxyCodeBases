#!/usr/bin/perl
use strict;
use warnings;
use threads;

my %file;
while (<>) {
  chomp;
  my $a = $_;
  s/\w+\///;
  if (open $file{$_}, "-|", "gzip -dc $a") {
    print STDERR "open $a was succesful!\n";
  } else {
    print STDERR "open $a failed!\n";
  }
}
print STDERR "Read file list complete!\n";

my %thread;
foreach (sort keys %file) {
  if (/sam.gz$/) {
    $thread{$_} = threads->new(\&countsam, $file{$_});
  } elsif (/fq.gz$/) {
    $thread{$_} = threads->new(\&countfq, $file{$_});
  }
}

my %count;
foreach (sort keys %thread) {
  $count{$_} = $thread{$_}->join;
  print STDERR "Count $_ complete!\n";
}

foreach (sort keys %count) {
  print $_, "\t", join( ',',@{$count{$_}} ), "\n";
}

sub countfq {
  my $in = shift;
  my ($lc,$f)=(0,0);
  while (<$in>) {
    ++$lc;
    if (/^\@[^ ]+ [12]:([YN]):/) {
    	++$f if $1 eq 'N';
    }
  }
  close $in;
  return [$lc/4,$f];
}

sub countsam {
  my $in = shift;
  my $lc;
  while (<$in>) {
    next if /^@/;
    ++$lc;
  }
  close $in;
  return [$lc/2];
}
