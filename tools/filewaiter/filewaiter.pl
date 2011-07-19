#!/bin/env perl
#use lib "/ifs1/ST_ASMB/USER/huxuesong/public/lib";
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.1;
our $opts='t:z';
our($opt_t, $opt_z);

#our $desc='';
our $help=<<EOH;
\t-t sleep seconds (3600)
\t-z file size can be zero
EOH
our $ARG_DESC='sampe_files';

ShowHelp();
$opt_t=3600 if ! $opt_t;
$opt_z=0 if ! $opt_z;

print STDERR "Waits changing of [@ARGV] for [$opt_t]s, Zero File[$opt_z]\n";

my $start_time = [gettimeofday];
#BEGIN
my %FileSizes0;
for my $name (@ARGV) {
    my $size=-1;
    $size = (-s $name) if (-f $name);
    $FileSizes0{$name} = $size;
    print "[$name]: $size\n";
}
my %FileSizes;
while (1) {
    sleep $opt_t;
    for my $name (@ARGV) {
        my $size=-1;
        $size = (-s $name) if (-f $name);
        $FileSizes{$name} = $size;
        print "[$name]: $size\n";
    }
    ;
}
#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__

