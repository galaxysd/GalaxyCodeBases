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

while (<>) {
    next if /^@\w\w\t\w\w:/;
    chomp;
    my @read1=split /\t/;
    <> or last;
    my @read2=split /\t/;
    die '[x]Not PE sam file.\n' if $read1[0] ne $read2[0];
    my ($QNAME,$FLAG,$RNAME,$POS,$MAPQ,$CIAGR,$MRNM,$MPOS,$ISIZE,$SEQ,$QUAL,$OPT)=@read1;
}

#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__

