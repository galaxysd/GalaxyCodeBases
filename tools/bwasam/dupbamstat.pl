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
        open SAM,'-|',"$SAMTOOLS view $bamfile" or (warn "[!]Error opening $bamfile: $!\n" and next);
        my ($READLEN,$i)=(0,0);
        while (defined($_=<SAM>) and $i<2048) {
            my @read1=split /\t/;
            next unless $read1[5] =~ /^(\d+)M$/;
            $READLEN = $1 if $READLEN < $1;
            ++$i;
        }
        close SAM;
        print STDERR "Read [$bamfile] $READLEN\n";

        open SAM,'-|',"$SAMTOOLS view $bamfile" or (warn "[!]Error opening $bamfile: $!\n" and next);
        while (<SAM>) {
                #next if /^(#|@)/;
                #chomp;
                my @read1=split /\t/;
                next unless $read1[1] & 3;  # paired + mapped in a proper pair
                next if $read1[1] >= 256;   # not primary || QC failure || optical or PCR duplicate
                next unless $read1[5] =~ /^(\d+)M$/;
                next unless $1 == $READLEN;
                next if $read1[11] eq 'XT:A:R'; # Type: Unique/Repeat/N/Mate-sw, N not found.
                #print "$_";
        }
        close SAM;
        #ddx \%NFO;
}

__END__

