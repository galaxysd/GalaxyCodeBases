#!/bin/env perl
=pod
Author: Hu Xuesong @ BGI <huxuesong@genomics.org.cn>
Version: 0.1.0 @ 20110803
=cut
#use lib "/ifs1/ST_ASMB/USER/huxuesong/public/lib";
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;
use Data::Dump qw(dump ddx);

$main::VERSION=0.1.0;
our $opts='o:b';
our($opt_o, $opt_b);

#our $desc='';
our $help=<<EOH;
\t-o output prefix (./matrixsum).{count,ratio}.matrix
\t-b No pause for batch runs
EOH
our $ARG_DESC='matrix_count_file';

ShowHelp();
$opt_o='./matrixsum' if ! $opt_o;
my $input=shift @ARGV;
die "[!]Error @ [$input]\n" unless (-s $input);

print STDERR "From [$input] to [$opt_o]\n";
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <STDIN>;}

my $start_time = [gettimeofday];
#BEGIN
open IN,'<',$input or die "Error: $!\n";
my $READLEN=0;
my $Qcount=41;
my ($TotalReads,$TotalBase,$MisBase,%BaseCountTypeRef)=(0,0,0);
my ($mapBase,$mapReads)=(0,0);
my $type='N/A';
my %Stat;   # $Stat{Ref}{Cycle}{Read}{Quality}
my %MismatchBYQ;   # $MismatchBYQ{Read:A,T,C,G,All}{$Q}->[mismatch,match]
my @BQHeader;
while (<IN>) {
    if (/^#Input \[(\w+)\] file of mapped Reads: (\d+) , mapped Bases (\d+) \(no base stat for sam files\)$/) {
        $type=$1 if $type ne 'sam';
        $mapReads += $2;
        $mapBase += $3 if $type ne 'sam';
    }
    if (/^#Total statistical Bases: (\d+) , Reads: (\d+) of ReadLength (\d+)$/) {
        print " >$input|   Reads:$2   ReadLen:$3   Bases:$1\n";
        $TotalReads += $2;
        $READLEN = $3 if $READLEN < $3;
    }
    if (/^#Ref\tCycle\t/) {
        #s/^#//;
        chomp;
        (undef,undef,@BQHeader)=split /\t/;
        pop @BQHeader if $BQHeader[-1] eq 'RowSum';
        for (@BQHeader) {
            /(\w)[^\d]?(\d+)$/ or die "[x]BQHeader wrong:[$_].\n";
            $_=[$1,$2];
        }
    }
    if (/^#Dimensions:.+?Quality_number (\d+)$/) {
        $Qcount = $1 if $Qcount<$1;
    }
    if (/^#Mismatch_base: (\d+)/) {
        $MisBase += $1;
    }
    next if /^#/;
    next if /^$/;
    chomp;
    my ($ref,$cycle,@BQ)=split /\t/;
    #print "$ref,$cycle,@BQ\n";
    next unless $ref =~ /^[ATCG]$/;
    #die "[$_]\n$ref,$cycle,[@BQ]\n$#BQ < $#BQHeader " if $#BQ < $#BQHeader;
    for my $key (@BQHeader) {
        my $value=shift @BQ;
        my ($kbase,$kQ)=@$key;
        $Stat{$ref}{$cycle}{$kbase}{$kQ}+=$value;
        if ($ref eq $kbase) {
            $MismatchBYQ{$kbase}{$kQ}->[1] +=$value;
            $MismatchBYQ{'_All'}{$kQ}->[1] +=$value;
        } else {
            $MismatchBYQ{$kbase}{$kQ}->[0] +=$value;
            $MismatchBYQ{'_All'}{$kQ}->[0] +=$value;
        }
        $BaseCountTypeRef{$ref}+=$value;
        $TotalBase+=$value;
        #print "{$ref}{$cycle}{$key}$value\n";
    }
}
close IN;
#ddx \%MismatchBYQ;
open OA,'>',$opt_o.'.err2mis' or die "Error: $!\n";
print OA "Read\tQ\tErrRate\tMismatchRate\terrbar\n";
for my $read (sort keys %MismatchBYQ) {
    for my $Q (sort {$a<=>$b} keys %{$MismatchBYQ{$read}}) {
        my ($mismatch,$match)=@{$MismatchBYQ{$read}{$Q}};
        next unless $match;
        print OA "$read\t$Q\t",10**(-$Q/10),"\t",$mismatch/($mismatch+$match),"\t",1/($mismatch+$match),"\n";
    }
}
close OA;
=pod
for my $ref (keys %Stat) {
    for my $cycle (keys %{$Stat{$ref}}) {
        for my $base (keys %{$Stat{$ref}{$cycle}}) {
            for my $Q (keys %{$Stat{$ref}{$cycle}{$base}}) {
                print "{$ref}{$cycle}{$base}{$Q} = $Stat{$ref}{$cycle}{$base}{$Q}\n";
            }
        }
    }
}
=cut

#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__

