#!/bin/env perl
#use lib "/ifs1/ST_ASMB/USER/huxuesong/public/lib";
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.1;
our $opts='r:o:b';
our($opt_o, $opt_r, $opt_b);

#our $desc='';
our $help=<<EOH;
\t-r ref fasta file (./ref/human.fa)
\t-o output prefix (./matrix).{raw,ratio}
\t-b No pause for batch runs
EOH
our $ARG_DESC='sampe_files';

ShowHelp();
$opt_r='./ref/human.fa' if ! $opt_r;
$opt_o='./matrix' if ! $opt_o;
die "[x]-r $opt_r not exists !\n" unless -f $opt_r;

print STDERR "From [@ARGV] with [$opt_r] to [$opt_o]\n";
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <STDIN>;}

my $start_time = [gettimeofday];
#BEGIN

warn "[!]Reading Reference Genome:\n";
my %Genome;
open GENOME,'<',$opt_r or die "Error: $!\n";
while (<GENOME>) {
	/^>(\S+)/ or next;
	my $seqname = $1;
    print STDERR " >$seqname ...";
	$/=">";
	my $genome=<GENOME>;
	chomp $genome;
	$genome=~s/\s//g;
	$/="\n";
    $Genome{$seqname}=$genome;
    print STDERR "\b\b\bdone !\n";
	$genome='';
}
close GENOME;
sub getBases($$$) {
    my ($chr,$start,$len)=@_;
    return substr $Genome{$chr},$start-1,$len;
}

my $READLEN=100;

my %Stat;   # $Stat{Ref}{Cycle}{Read}{Quality}
sub statRead($$$$$) {
    my ($ref,$isReverse,$read,$Qstr,$cyclestart)=@_;
    my $PEpos;
    for (my $i=0;$i<$READLEN;$i++) {
        my $refBase=substr $ref,$i,1;
        my $readBase=substr $read,$i,1;
        my $QstrSingle=substr $Qstr,$i,1;
        my $Qval=ord($QstrSingle)-33;
        if ($isReverse) {
            $PEpos=$cyclestart+$READLEN-1-$i;
        } else {
            $PEpos=$cyclestart+$i;
        }
        ++$Stat{$refBase}{$PEpos}{$readBase}{$Qval};
        
print "$isReverse {$refBase}{$PEpos}{$readBase}{$Qval} ",($refBase eq $readBase)?'=':'x',"\n";
    }
}

while (<>) {
    next if /^@\w\w\t\w\w:/;
    chomp;
    my @read1=split /\t/;
    chomp($_=<>) or last;
    my @read2=split /\t/;
print join("\t",@read1),"\n-",join("\t",@read2),"\n";
    die "[x]Not PE sam file.\n" if $read1[0] ne $read2[0];
    next unless $read1[1] & 3;
    next if $read1[1] >= 256;
    next unless $read2[1] & 3;
    next if $read2[1] >= 256;
    next unless $read1[5] =~ /^(\d+)M$/;
    next unless $1 == $READLEN;
    next unless $read2[5] =~ /^(\d+)M$/;
    next unless $1 == $READLEN;
    next unless $read1[6] eq '=';
    next unless $read2[6] eq '=';
    next if $read1[11] eq 'XT:A:R';
    next if $read2[11] eq 'XT:A:R';
    my $ref1=uc getBases($read1[2],$read1[3],$READLEN) or print join("\t",@read1),"\n";
    my $ref2=uc getBases($read2[2],$read2[3],$READLEN) or print join("\t",@read2),"\n";
    #my ($QNAME,$FLAG,$RNAME,$POS,$MAPQ,$CIAGR,$MRNM,$MPOS,$ISIZE,$SEQ,$QUAL,$OPT)=@read1;
    #       0      1    2       3   4       5   6       7     8     9    10    11
    statRead($ref1,$read1[1] & 16,$read1[9],$read1[10],1);
    statRead($ref2,$read2[1] & 16,$read2[9],$read2[10],101);
}

#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__

