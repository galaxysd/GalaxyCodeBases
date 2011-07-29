#!/bin/env perl
#use lib "/ifs1/ST_ASMB/USER/huxuesong/public/lib";
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.1;
our $opts='o:r:w:m:n:b';
our($opt_o, $opt_r, $opt_w, $opt_m, $opt_n, $opt_b);

#our $desc='';
our $help=<<EOH;
\t-w window sizes (100,200,500,1000)
\t-r ref fasta file (./ref/human.fa)
\t-o output prefix (./gcdepdat).{mcount,mratio}
\t-m min win non-N-base count (30)
\t-n max N ratio in one win (0.75)
\t-b No pause for batch runs
EOH
our $ARG_DESC='coverage_fa_files{,.gz,.bz2}';

ShowHelp();
$opt_w='100,200,500,1000' if ! $opt_w;
$opt_r='./ref/human.fa' if ! $opt_r;
$opt_m=30 if ! $opt_m;
$opt_n=0.75 if ! $opt_n;
$opt_o='./gcdepdat' if ! $opt_o;
die "[x]No input files found !\n" unless @ARGV;
die "[!]Max 252 files supported.\n" if @ARGV>252;

my @wins=grep {$_>=50} map {int $_} split /,/,$opt_w;
die "[x]Window Size must >= 50.\n" unless @wins;

print STDERR "From [@ARGV] with [$opt_r][$opt_m][$opt_n] to [$opt_o] of [",join(',',@wins),"]\n";
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

my (%Genome,%EffChrLen);
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
    my $n=($genome=~s/[^ATCG]/A/ig);
    $EffChrLen{$seqname}=length($Genome{$seqname})-$n;
    print STDERR "\b\b\b   \t",length $Genome{$seqname},".\n";
	$genome='';
}
close GENOME;

my (%Result,%Stat,%CVG,%Depth,%DepCnt);

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
	            $DepDatChr[$i]+=$1 if ($1!=65536);
	            ++$i;
	        }
	        $genome='';
	        last;
	    }
	    print STDERR '.';
	}
	for (@DepDatChr) {
	    $_=65535 if $_>65535;
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
	        my $size=$win-$n;
	        next if $size<$opt_m or ($n/$win)>$opt_n;
	        $gc=int($gc/$size);
	        my $sum=0;
	        $sum+=$DepDatChr[$_] for ($start..$start+$win);
	        my $value=$sum/$size;
	        push @{$Result{$win}{$gc}},$value;
	        $Stat{$win}[$gc][0] += $value;
	        ++$Stat{$win}[$gc][1];
	        $start += $win;
	    }
	}
	for (@DepDatChr) {
	    if ($_) {
	        ++$CVG{$SeqName};
	        ++$CVG{'__ALL__'};
	        $Depth{$SeqName}+=$_;
	        $Depth{'__ALL__'}+=$_;
	        ++$DepCnt{$_}{$SeqName};
	        ++$DepCnt{$_}{'__ALL__'};
	    }
	}
}
close $_ for (@FH,$firstFH);
@FH=();
my @Chr=sort keys %Genome;

open ODEP,'>',$opt_o.'_stat.tsv' or die "Error: $!\n";
print ODEP "#ChrID\tDepth\tCovered\tCVGratio\tEffChrLen\tNzone\tChrLen\n";
for my $chr ('__ALL__',@Chr) {
    print ODEP join("\t",$chr,$Depth{$chr}/$EffChrLen{$chr},$CVG{$chr},$CVG{$chr}/$EffChrLen{$chr},$EffChrLen{$chr},length($Genome{$chr})-$EffChrLen{$chr},length($Genome{$chr})),"\n";
}
close ODEP;

open OHST,'>',$opt_o.'_hist.tsv' or die "Error: $!\n";
print OHST "#Depth\t__ALL__";
print OHST "\t$_" for @Chr;
print OHST "\n";
for my $dep (sort {$a<=>$b} keys %DepCnt) {
    print OHST $dep;
    for my $chr ('__ALL__',@Chr) {
        my $v=0;
        $v=$DepCnt{$dep}{$chr} if exists $DepCnt{$dep}{$chr};
        print OHST "\t",$v;
    }
    print OHST "\n";
}
close OHST;

for my $win (@wins) {
    open O,'>',$opt_o."_gcdep.$win.tsv" or die "Error: $!\n";
    print O "#WinSize=$win\tWinCount=$Stat{$win}[1]
#GC%\tRefCnt\tMean\tSmall\tQ1\tMiddle\tQ3\tBig\tmin\tmax\trefcntcal\n";
    for my $gc (sort {$a<=>$b} keys %{$Result{$win}}) {
        print O $gc,'.5',"\t",$Stat{$win}[$gc][1],"\t",$Stat{$win}[$gc][0]/$Stat{$win}[$gc][1],"\t";
        ###
        print O "\n";
    }
}


#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__

