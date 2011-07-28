#!/bin/env perl
#use lib "/ifs1/ST_ASMB/USER/huxuesong/public/lib";
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.1;
our $opts='o:r:w:b';
our($opt_o, $opt_r, $opt_w, $opt_b);

#our $desc='';
our $help=<<EOH;
\t-w window sizes (100,200,500,1000)
\t-r ref fasta file (./ref/human.fa)
\t-o output prefix (./gcdepdat).{mcount,mratio}
\t-b No pause for batch runs
EOH
our $ARG_DESC='coverage_fa_files{,.gz,.bz2}';

ShowHelp();
$opt_w='100,200,500,1000' if ! $opt_w;
$opt_r='./ref/human.fa' if ! $opt_r;
$opt_o='./gcdepdat' if ! $opt_o;
die "[x]No input files found !\n" unless @ARGV;
die "[!]Max 252 files supported.\n" if @ARGV>252;

my @wins=grep {$_>=50} map {int $_} split /,/,$opt_w;
die "[x]Window Size must >= 50.\n" unless @wins;

print STDERR "From [@ARGV] with [$opt_r] to [$opt_o] of [",join(',',@wins),"]\n";
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <STDIN>;}

my $start_time = [gettimeofday];
#BEGIN
my @FH;
while($_=shift @ARGV) {
    my $infile;
    if (/.bz2$/) {
	    open( $infile,"-|","bzip2 -dc $_") or die "Error opening $_: $!\n";
    } elsif (/.gz$/) {
     	open( $infile,"-|","gzip -dc $_") or die "Error opening $_: $!\n";
    } else {open( $infile,"<",$_) or die "Error opening $_: $!\n";}
    push @FH,$infile;
}
warn '[!]depth files opened: ',scalar @FH,"\n[!]Reading Reference Genome:\n";

my %Genome;
open GENOME,'<',$opt_r or die "Error: $!\n";
while (<GENOME>) {
    s/^>//;
	/^(\S+)/ or next;
	my $seqname = $1;
    print STDERR " >$seqname ...";
	$/=">";
	my $genome=<GENOME>;
	chomp $genome;
	$genome=~s/\s//g;
	$/="\n";
    $Genome{$seqname}=$genome;
    print STDERR "\b\b\b   \t",length $Genome{$seqname},".\n";
	$genome='';
}
close GENOME;

my %Result;

warn "[!]Reading Depth Files:\n";
my $firstFH=shift @FH;
while(<$firstFH>) {
    s/^>//;
	/^(\S+)/ or next;
	my $SeqName = $1;
    print STDERR " \@$SeqName ...";
	$/=">";
	my $genome=<$firstFH>;
	$/="\n";
	#print STDERR length $genome,"|\t";
    my @DepDatChr=();
    while($genome=~/(\d+)/g) {
        my $v=$1;
        $v=0 if $v==65536;
        push @DepDatChr,$v;
    }
    print STDERR "\b\b\b   \t",scalar @DepDatChr,': .';
	$genome='';
	for my $fh (@FH) {
        while(<$fh>) {
            s/^>//;
	        /^(\S+)/ or next;
	        die "[x]Depth file in different order ($SeqName ne $1) !\n" if $SeqName ne $1;
	        $/=">";
	        $genome=<$fh>;
	        $/="\n";
	        my $i=0;
	        while($genome=~/(\d+)/gc) {
	            $DepDatChr[$i]+=$1 if ($1!=65536) and ($DepDatChr[$i]+$1<65535);
	            ++$i;
	        }
	        $genome='';
	        last;
	    }
	    print STDERR '.';
	}
	print STDERR "\n";
	# one chr done
	for my $win (@wins) {
	    my $chrlenOK=length($Genome{$SeqName})-$win;
	    my $start=0;
	    while($start<$chrlenOK) {   # the last win is canceled ...
	        my $seq=substr $Genome{$SeqName},$start,$win;
	        my $gc=($seq=~s/[GC]/A/ig);
	        my $n=($seq=~s/[^ATCG]/A/ig);
	        my $size-$win-$n;
	        $gc=int($gc/$size);
	        my $sum=0;
	        $sum+=$DepDatChr[$_] for ($start..$start+$win);
	        push @{$Result{$win}->[$gc]},$sum/$size;
	        $start += $win;
	    }
	}
}



close $_ for (@FH,$firstFH);
#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__

