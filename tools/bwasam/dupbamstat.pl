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
        print "perl $0 <outfile> <sorted bam files>\n";        # soaps.nfo can store the file size of soap. Too
        exit;
}

my $statout = shift @ARGV;
my ($DupSE,$DupPE,$ReadsStat)=(0,0,0);
my $FilesStr='[' . join('], [',@ARGV) . "]";

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

        my ($lastL,$lastR,$readL,$readR,$isDupSE,$isDupPE)=(0,0,0,0,0,0);
        open SAM,'-|',"$SAMTOOLS view $bamfile" or (warn "[!]Error opening $bamfile: $!\n" and next);
        while (<SAM>) {
                #next if /^(#|@)/;
                #chomp;
                my @read1=split /\t/;
                next if $read1[2] eq '*';   # unmap
                next unless $read1[1] & 3;  # paired + mapped in a proper pair
                next if $read1[1] >= 256;   # not primary || QC failure || optical or PCR duplicate
                next unless $read1[5] =~ /^(\d+)M$/;
                next unless $1 == $READLEN;
                next if $read1[11] eq 'XT:A:R'; # Type: Unique/Repeat/N/Mate-sw, N not found.
                #print "$_";
                if ($lastL!=0) {
                    ($readL,$readR)=@read1[3,7];
                    if ($readL == $lastL) {
                        $isDupSE=1;
                        if ($readR == $lastR) {
                            $isDupPE=1;
                        }
                    } else {
                        ($lastL,$lastR)=($readL,$readR);
                        $DupSE += $isDupSE;
                        $DupPE += $isDupPE;
                        ($isDupSE,$isDupPE)=(0,0);
                    }
                    ++$ReadsStat;
                    next;
                } else {
                    ($lastL,$lastR)=@read1[3,7];
                    next;
                }
=pod
                         ┌────┬───────┬──────────────────────────────────────────────────────────┐
                         │Col │ Field │                       Description                        │
                         ├────┼───────┼──────────────────────────────────────────────────────────┤
                         │ 1  │ QNAME │ Query (pair) NAME                                        │
                         │ 2  │ FLAG  │ bitwise FLAG                                             │
                         │ 3  │ RNAME │ Reference sequence NAME                                  │
                         │ 4  │ POS   │ 1-based leftmost POSition/coordinate of clipped sequence │
                         │ 5  │ MAPQ  │ MAPping Quality (Phred-scaled)                           │
                         │ 6  │ CIAGR │ extended CIGAR string                                    │
                         │ 7  │ MRNM  │ Mate Reference sequence NaMe (`=' if same as RNAME)      │
                         │ 8  │ MPOS  │ 1-based Mate POSistion                                   │
                         │ 9  │ ISIZE │ Inferred insert SIZE                                     │
                         │10  │ SEQ   │ query SEQuence on the same strand as the reference       │
                         │11  │ QUAL  │ query QUALity (ASCII-33 gives the Phred base quality)    │
                         │12  │ OPT   │ variable OPTional fields in the format TAG:VTYPE:VALUE   │
                         └────┴───────┴──────────────────────────────────────────────────────────┘
=cut
        }
        close SAM;
        #ddx \%NFO;
}
open OA,'>',$statout or die "Error: $!\n";
print OA "Input Files: $FilesStr\nTotal Reads: $ReadsStat/2 = ",$ReadsStat/2,"\nSE Dup: ",$DupSE/2,"\nPE dup: ",$DupPE/2,"\n",'-'x80,"\nPE dup ratio: ",$DupPE/$ReadsStat,"\nSE dup ratio: ",$DupSE/$ReadsStat,"\n";
close OA;

__END__

