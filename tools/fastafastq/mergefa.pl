#!/bin/env perl
#use lib "/ifs1/ST_ASMB/USER/huxuesong/public/lib";
use lib '/export/data0/gentoo/develop/toGit/perl/perlib/etc';
use strict;
use warnings;
use Time::HiRes qw ( gettimeofday tv_interval );
use Galaxy::ShowHelp;

$main::VERSION=0.0.1;
our $opts='o:cb';
our($opt_o, $opt_c, $opt_b);

#our $desc='';
our $help=<<EOH;
\t-o output file
\t-c keep only ChrID ~ /^chr\\d+\$/
\t-b No pause for batch runs
EOH
our $ARG_DESC='fa_files{,.gz,.bz2}';

ShowHelp();
die "[x]No input files found !\n" unless @ARGV;
die "[!]Max 252 files supported.\n" if @ARGV>252;

die "[x]Need output file !\n" unless $opt_o;

print STDERR "From [@ARGV] to [$opt_o]\n";
print STDERR "Keep only ChrID =~ /^chr\\d+\$/\n" if $opt_c;
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
warn '[!]Files opened: ',scalar @FH,"\n[!]Reading:\n";

my $OutFile;
open $OutFile,'>',$opt_o  or die "Error opening $opt_o: $!\n";
print $OutFile ">Seq\n";

sub rwfa($$) {
	my ($ifh,$ofh)=@_;
	while (<$ifh>) {
		s/^>//;
		/^(\S+)/ or next;
		my $seqname = $1;
		if ($opt_c) {
			next if $seqname !~ /^chr\d+$/;
		}
		print STDERR ">$seqname\n";
		$/=">";
		my $genome=<$ifh>;
		chomp $genome;
		$genome=~s/[^ATCGatcg]//g;
		$/="\n";
		print $OutFile $genome;
		$genome='';
	}
	close $ifh;
}

&rwfa($_,$OutFile) for @FH;
print $OutFile "\n";
close $OutFile;



#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
__END__
my (%Genome,%EffChrLen);
sub readfa($) {
	my $fh=$_[0];
	while (<$fh>) {
		s/^>//;
		/^(\S+)/ or next;
		my $seqname = $1;
		#print STDERR " >$seqname ...";
		$/=">";
		my $genome=<$fh>;
		chomp $genome;
		#$genome=~s/\s//g;
		$genome=~s/[^ATCG]//g;
		$/="\n";
		$Genome{$seqname}=$genome;
		#my $n=($genome=~s/[^ATCG]/A/ig);
		#$EffChrLen{$seqname}=length($Genome{$seqname})-$n;
		#print STDERR "\b\b\b   \t",length $Genome{$seqname},".\n";
		$genome='';
	}
	close $fh;
}

