#!/bin/env perl
use strict;
use warnings;
#use Data::Dump qw(dump ddx);

chomp(my $user=`id -nru`);
my ($SAMTOOLS);
if ($user eq 'galaxy') {
    $SAMTOOLS='samtools';
} else {
    $SAMTOOLS='/ifs1/ST_ASMB/USER/huxuesong/group/illumia/_reads/samtools';
}

unless (@ARGV){
        print "perl $0 <outprefix(.stat & .ins)> <sorted bam files>\n";        # soaps.nfo can store the file size of soap. Too
        exit;
}

my $statout = shift @ARGV;

while(my $bamfile=shift @ARGV) {
    next unless -f $bamfile;

        print STDERR "Read [$bamfile]\n";

        open SAM,'-|',"$SAMTOOLS view $bamfile" or (warn "[!]Error opening $bamfile: $!\n" and next);
        while (<SAM>) {
                next if /^(#|$)/;
                chomp;
        }
        close SAM;
        #ddx \%NFO;
}

__END__

