#!/bin/env perl
use strict;
use warnings;
use lib '/share/raid010/resequencing/soft/lib';
use Time::HiRes qw ( gettimeofday tv_interval );
use File::Basename;
use Galaxy::ShowHelp;

$main::VERSION=0.0.1;

our $opts='i:r:g:p:n:s:o:vb';
our($opt_i, $opt_r, $opt_g, $opt_p, $opt_n, $opt_s, $opt_v, $opt_b, $opt_o);

our $desc='';
our $help=<<EOH;
\t-i Genome Info [ChrID\\tLen]
\t-r RegEx for ChrID (chromosome_\\d+)
\t-g GFF file for CDS length
\t-p SNP file
\t-n Indel file
\t-s SV file
\t-o Output file
\t-v show verbose info to STDOUT
\t-b No pause for batch runs
EOH

ShowHelp();

$opt_r = 'chromosome_\d+' unless $opt_r;

print STDERR "From [$opt_i][$opt_g][$opt_p][$opt_n][$opt_s] to [$opt_o] with to [$opt_r]\n";
unless ($opt_b) {print STDERR 'press [Enter] to continue...'; <>;}

my $start_time = [gettimeofday];
#BEGIN
open( IN,'<',$opt_i) or die "[x]Error: $!\n";
my $GenomeLen;
while (<IN>) {
	chomp;
	my ($chrid,$len)=split /\t/;
	next unless /$opt_r/;
	$GenomeLen += $len;
	print "$chrid,$len\t$GenomeLen\n" if $opt_v;
}
close IN;

#END
my $stop_time = [gettimeofday];

print STDERR "\nTime Elapsed:\t",tv_interval( $start_time, $stop_time )," second(s).\n";
#print "\nPlease use the following command to batch qsub:\033[32;1m
#find ${opt_o}/ -name '*.sh' | while read ll; do qsub -l vf=2G -cwd \$ll; done\n\033[0;0m\n";
__END__
