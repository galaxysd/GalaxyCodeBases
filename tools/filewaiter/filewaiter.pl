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

print STDERR "Waits changing of [@ARGV] for [$opt_t]s, Zero File[",$opt_z?'OK':'Wait',"]\n";

my $start_time = [gettimeofday];
#BEGIN
my (%FileSizes0,%Files);
for my $name (@ARGV) {
    my $size=-1;
    $size = (-s $name) if (-f $name);
    $FileSizes0{$name} = $size;
    ++$Files{$name};
    print "[$name]: $size\n";
}
my %FileSizesC;
while (keys %Files) {
    sleep $opt_t;
    for my $name (keys %Files) {
        my $size=-1;
        chomp(my $date=`date`);
        if (-f $name) {
            $size = (-s $name);
        } else {
            print "$date [!]$name. NULL\n";
            next;
        }
        unless ($opt_z) {
            if ($size==0) {
                print "$date [d]$name. $FileSizes0{$name} -> $size\n";
                $FileSizes0{$name} = $size;
                next;
            }
        }
        if ($size != $FileSizes0{$name}) {
            print "$date [d]$name. $FileSizes0{$name} -> $size\n";
            $FileSizes0{$name} = $size;
        } else {
            ++$FileSizesC{$name};
            if ($FileSizesC{$name} > 2) {
                delete $Files{$name};
                print "$date [O]$name. Stablie @ $size\n\n";
            } else {
                print "$date [s]$name. $FileSizes0{$name} == $size\n";
            }
        }
    }
}
#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__

