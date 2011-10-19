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
        print "perl $0 <maxPEpairs(0=Inf.)> <outfile> <sorted bam files>\n";        # soaps.nfo can store the file size of soap. Too
        exit;
}

my ($STAGELENGTH,$stageborder)=(2000,2000);

my $maxpairs = shift @ARGV;
my $statout = shift @ARGV;
my ($DupSE,$DupPE,$ReadsStat,%SEDup,%PEDup)=(0,0,0);
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

        my ($lastID,$ReadID,$lastL,$lastR,$readL,$readR,$isDupSE,$isDupPE)=('','',0,0,0,0,0,0);
        my ($MaxIns,$ins)=(0);
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
                    ($ReadID,$readL,$readR,$ins)=@read1[0,3,7,8];
                    $MaxIns = $ins if $MaxIns < $ins;
                    ++$SEDup{$readL}{$ReadID};
                    my ($reada,$readb) = sort {$a<=>$b} ($readL,$readR);
                    ++$PEDup{$readb}{$reada}{$ReadID};
                    if ($reada>$stageborder) {
                        for (sort {$a<=>$b} keys %SEDup) {
                            my $dupcnt=scalar keys $SEDup{$_}->{};
                            $DupSE += $dupcnt if $dupcnt>=2;
                        }
                        %SEDup=();
                        for my $Rb (sort {$a<=>$b} keys %PEDup) {
                            if ($Rb<=$stageborder) {
                                for my $Ra (sort {$a<=>$b} keys %{$PEDup{$Rb}}) {
                                    my $dupcnt=scalar keys $PEDup{$Rb}{$Ra}->{};
                                    $DupPE += $dupcnt if $dupcnt>=2;
                                }
                                delete $PEDup{$Rb};
                            }
                        }
                        $stageborder += $STAGELENGTH;
                    }
=pod
                    if ($readL == $lastL) {
                        $isDupSE=1;
                        ++$SEDup{$ReadID};
                        if ($readR == $lastR) {
                            $isDupPE=1;
                            ++$PEDup{$ReadID};
                        }
                    } else {
                        if ($isDupSE) {
                            ++$DupSE;
                            ++$PEDup{$lastID};
                        }
                        if ($isDupPE) {
                            ++$DupPE;
                            ++$SEDup{$lastID};
                        }
                        ($lastID,$lastL,$lastR)=($ReadID,$readL,$readR);
                        ($isDupSE,$isDupPE)=(0,0);
                    }
                    ++$ReadsStat;
                    next;
=cut                    
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

