#!/bin/env perl
#use lib "/ifs1/ST_ASMB/USER/huxuesong/public/lib";
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.1;
our $opts='o:b';
our($opt_o, $opt_b);

#our $desc='';
our $help=<<EOH;
\t-o output prefix (./allmatrix).{mcount,mratio}
\t-b No pause for batch runs
EOH
our $ARG_DESC='matrix_count_files';

ShowHelp();
$opt_o='./allmatrix' if ! $opt_o;

print STDERR "From [@ARGV] to [$opt_o]\n";
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <STDIN>;}

my $start_time = [gettimeofday];
#BEGIN
my $READLEN=0;
my ($TotalReads,$TotalBase,%BaseCountTypeRef)=(0);
my %Stat;   # $Stat{Ref}{Cycle}{Read-Quality}
my @BQHeader;
while (<>) {
    if (/^#Total statistical Bases: (\d+) , Reads: (\d+) of ReadLength (\d+)$/) {
        print " >$ARGV|   Reads:$2   ReadLen:$3   Bases:$1\n";
        $TotalReads += $2;
        $READLEN = $3 if $READLEN < $3;
    }
    if (/^#Ref\tCycle\t/) {
        #s/^#//;
        chomp;
        (undef,undef,@BQHeader)=split /\t/;
        pop @BQHeader if $BQHeader[-1] eq 'RowSum';
    }
    next if /^#/;
    next if /^$/;
    chomp;
    my ($ref,$cycle,@BQ)=split /\t/;
    #print "$ref,$cycle,@BQ\n";
    #die "[$_]\n$ref,$cycle,[@BQ]\n$#BQ < $#BQHeader " if $#BQ < $#BQHeader;
    for my $key (@BQHeader) {
        my $value=shift @BQ;
        $Stat{$ref}{$cycle}{$key}+=$value;
        $BaseCountTypeRef{$ref}+=$value;
        $TotalBase+=$value;
        #print "{$ref}{$cycle}{$key}$value\n";
    }
}
#print $TotalReads,"\t",$READLEN,"\n";
#print join("\t",@BQHeader),"\n";
open OA,'>',$opt_o.'.mcount' or die "Error: $!\n";
open OB,'>',$opt_o.'.mratio' or die "Error: $!\n";
my $tmp;
$tmp="#Total statistical Bases: $TotalBase , Reads: $TotalReads of ReadLength $READLEN
#Reference Base Ratio in reads: ";
my @BaseOrder=sort keys %BaseCountTypeRef;  # qw{A T C G};
for (@BaseOrder) {
    $tmp .= $_.' '. int(0.5+100*1000*$BaseCountTypeRef{$_}/$TotalBase)/1000 .' %;   ';
}

$tmp .= "\n#".join("\t",'Ref','Cycle',@BQHeader);
print OA $tmp;
print OB $tmp;
print OA "\tRowSum";
print OB "\n";
my ($count,$countsum);
for my $ref (@BaseOrder) {
    print OA "\n";
    for my $cycle (sort {$a<=>$b} keys %{$Stat{$ref}}) {
        $tmp="$ref\t$cycle\t";
        print OA $tmp; print OB $tmp;
        my (@Counts,@Rates)=();
        for my $bq (@BQHeader) {
            push @Counts,$Stat{$ref}{$cycle}{$bq};
        }
        #print "[",join("|",@Counts),"\n";
        $countsum=0;
        $countsum += $_ for @Counts;
        push @Rates,$_/$countsum for @Counts;
        print OA join("\t",@Counts,$countsum),"\n";
        print OB join("\t",@Rates),"\n";
    }
}
close OA;
close OB;

#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__

