#!/bin/env perl
use lib '/nas/RD_09C/resequencing/soft/lib';
use lib 'E:/BGI/toGit/perlib/etc';
use strict;
use warnings;
use File::Basename;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.1;


our $opts='i:fb';
our($opt_i, $opt_f, $opt_b);

our $help=<<EOH;
\t-i Input FASTA file
\t-f Also Split to [fabychr/]
\t-b No pause for batch runs
EOH

ShowHelp();

die "[x]Must specify FASTA file !\n" if ! $opt_i;
die "[x]-i $opt_i not exists !\n" unless -f $opt_i;

my $FASTAwidth=80;

my ($filename, $outpath) = fileparse($opt_i);
unless (-w $outpath) {
	warn "[!]Cannot write to [$outpath]. Change to [./].\n";
	$outpath='./';
}

my $cutfa;
unless ($opt_f) {
	$cutfa=", Stat. Only.";
} else {
	$opt_f = 'fabychr';
	$cutfa=" to [$outpath$opt_f/]";
}

print STDERR "From [$opt_i]$cutfa\n";
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <STDIN>;}

my $start_time = [gettimeofday];
if ($opt_f) {
	mkdir $outpath.$opt_f,0755;
	die "[x]Error creating [$outpath$opt_f]. $!\n" if $! && $! ne 'File exists';
}

#BEGIN
my $infile;
if ($opt_i =~ /.bz2$/) {
	open( $infile,"-|","bzip2 -dc $opt_i") or die "Error: $!\n";
} elsif ($opt_i =~ /.gz$/) {
 	open( $infile,"-|","gzip -dc $opt_i") or die "Error: $!\n";
} else {open( $infile,"<",$opt_i) or die "Error: $!\n";}

{
	local $/=">";
	$_=<$infile>;
	die "[x]Not a FASTA file !\n" unless /^\s*>/;
}

open CHRO,'>',$outpath.'chrorder' or die "Error: $!\n";
open STAT,'>',$outpath.'chr.nfo' or die "Error: $!\n";
print STAT '#',join("\t",qw/ChrID Len EffLen GCratio N n Xx lc GC ChrDesc/),"\n";
while (<$infile>) {
	chomp;
	my $Head;
	my ($id,$desc)=split / /,$_,2;
	if ($desc && $desc !~ /^\s*$/) {
		$desc=~s/\t/_/g;
		$Head="$id $desc";
	} else { $desc='.';$Head=$id; }
	$/=">";
	my $seq=<$infile>;
	chomp $seq;
	$seq =~ s/\s//g;
	$/="\n";
	print STDERR ">$id,\t[$desc]:\n";
	print CHRO "$id\n";
	my $len=length($seq);
	if ($opt_f) {
		open O,'>',"$outpath$opt_f/$id.fa" or die "Error: $!\n";
		print O ">$Head\n";
		for (my $i=0; $i<$len; $i+=$FASTAwidth) {
			print O substr($seq,$i,$FASTAwidth)."\n";
		}
		close O;
	}
	my $N = $seq=~tr/N//;	# stand for gap or low quality
	my $n = $seq=~tr/n//;
	my $lc = $seq=~tr/agct/AGCT/;	# may stand for masked region
	my $GC = $seq=~tr/GCgc//;
	my $X = $seq=~tr/Xx//;	# stand for masked region
	my $Efflen=$len-$N-$n;
	my $GCratio=sprintf('%.3f',int(0.5+1000*$GC/$len)/1000);
	print STAT join("\t",$id,$len,$Efflen,$GCratio,$N,$n,$X,$lc,$GC,$desc),"\n";
	print STDERR " Len:$len, Effictive_Len:$Efflen, GC_Ratio:$GCratio\n";
}
close CHRO;
close STAT;

#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__
